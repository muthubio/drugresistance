#!/usr/bin/env python3
"""
pipeline_core.py
Core library for: Trim Galore -> BWA (paired + singles) -> merge -> fixmate -> sort -> markdup -> index
Optional: GATK HaplotypeCaller (+ bcftools normalization)
NEW:
  - --mask-bed: exclude problematic regions before annotation (like TB-Profiler)
  - --dr-catalog: annotate drug resistance using a JSON catalogue (e.g., tbdb.dr.json)
  - cDNA/HGVS support via snpEff (match 'RvXXXX' + 'c.###N>M')
  - Outputs: <sample>.dr_mutations.tsv, <sample>.dr_summary.tsv, <sample>.tbprofiler_like.json
"""
from __future__ import annotations

import re, shlex, shutil, subprocess as sp, tempfile, time, json, csv
from dataclasses import dataclass
from pathlib import Path
from typing import Optional, Tuple, Dict, List, Any, DefaultDict
from uuid import uuid4
from contextlib import contextmanager
from collections import defaultdict

# ---------------- logging & command runner (TB-Profiler style) ----------------
_LEVELS = {"DEBUG": 10, "INFO": 20, "WARNING": 30, "ERROR": 40}
_LOG_LEVEL = "INFO"
_LOG_FH = None  # file handle if --log-file is set

def _ts() -> str: return time.strftime("[%H:%M:%S]")

def set_log_level(level: str):
    global _LOG_LEVEL
    lvl = (level or "INFO").upper()
    if lvl not in _LEVELS: lvl = "INFO"
    _LOG_LEVEL = lvl

def set_log_file(path: str | None):
    global _LOG_FH
    if _LOG_FH:
        try: _LOG_FH.close()
        except: pass
    _LOG_FH = open(path, "a") if path else None

def _emit(level: str, msg: str):
    line = f"{_ts()} {level:<8} {msg}"
    print(line, flush=True)
    if _LOG_FH:
        _LOG_FH.write(line + "\n"); _LOG_FH.flush()

def DEBUG(msg: str):
    if _LEVELS[_LOG_LEVEL] <= _LEVELS["DEBUG"]:
        _emit("DEBUG", msg)

def INFO(msg: str):
    if _LEVELS[_LOG_LEVEL] <= _LEVELS["INFO"]:
        _emit("INFO", msg)

def WARNING(msg: str): _emit("WARNING", msg)
def ERROR(msg: str):   _emit("ERROR", msg)

def close_logs():
    global _LOG_FH
    if _LOG_FH:
        try: _LOG_FH.close()
        finally: _LOG_FH = None

def run_cmd(cmd, shell: bool = False, capture_stdout: bool = False, workdir: str | None = None):
    if isinstance(cmd, list): cmd_str = " ".join(shlex.quote(str(c)) for c in cmd)
    else: cmd_str = str(cmd)
    DEBUG(f"Running command: {cmd_str: <100} utils.py:480")
    proc = sp.run(
        cmd, shell=shell, cwd=workdir,
        stdout=sp.PIPE if capture_stdout else None,
        stderr=sp.PIPE, text=True
    )
    if proc.returncode != 0:
        ERROR("Command Failed:\n" + cmd_str + f"\nstderr:\n{proc.stderr or ''}")
        raise SystemExit(1)
    return proc.stdout if capture_stdout else None

@contextmanager
def step(title: str, src: str = ""):
    tail = f"  {src}" if src else ""
    INFO(f"{title}{tail}")
    t0 = time.time()
    try:
        yield
        DEBUG(f"Finished: {title} in {time.time()-t0:.1f}s{tail}")
    except Exception as e:
        ERROR(f"Failed: {title}{tail}  ({e})")
        raise

# ---------------- small helpers ----------------
def which_or_die(name: str) -> str:
    p = shutil.which(name)
    if p is None:
        raise SystemExit(f"ERROR: required program '{name}' not found in PATH")
    return p

def require_files(*paths: str) -> None:
    for p in paths:
        if not Path(p).exists():
            raise SystemExit(f"ERROR: required file not found: {p}")

def gz_cat_to(out_path: Path, *in_paths: Optional[Path]) -> None:
    out_path = Path(out_path)
    with open(out_path, "wb") as w:
        for ip in in_paths:
            if not ip: continue
            p = Path(ip)
            if not p.exists() or p.stat().st_size == 0:
                continue
            with open(p, "rb") as r:
                shutil.copyfileobj(r, w)

# ---------------- config models ----------------
@dataclass
class PipelineConfig:
    # IO
    r1: str
    r2: str
    ref: str
    sample: str
    out_bam: Optional[str] = None
    tmpdir: Optional[str]  = None
    sample_ploidy: int = 1

    # mapping/trimming
    length: int = 36
    threads: int = 1
    mem: str = "768M"
    fastqc: bool = False
    skip_singles: bool = False

    # calling
    call_variants: bool = False
    emit_gvcf: bool = False
    intervals: Optional[str] = None
    normalize: bool = False
    gatk_extra: Optional[str] = None

    # DR annotation
    dr_catalog: Optional[str] = None  # JSON file (e.g., tbdb.dr.json)
    min_alt_af: float = 0.0
    min_dp: int = 0
    mask_bed: Optional[str] = None    # BED to exclude (like TB-Profiler)

    # snpEff (for cDNA/HGVS matching)
    snpeff_db: Optional[str] = "Mycobacterium_tuberculosis_h37rv"

@dataclass
class TrimOutputs:
    out_dir: Path
    val1: Path
    val2: Path
    unp1: Optional[Path]
    unp2: Optional[Path]

# ---------------- core steps ----------------
def prepare_tempdir(cfg: PipelineConfig) -> Tuple[Path, bool]:
    if cfg.tmpdir:
        tmp = Path(cfg.tmpdir); tmp.mkdir(parents=True, exist_ok=True)
        return tmp, True
    tmp_ctx = tempfile.TemporaryDirectory()
    return Path(tmp_ctx.name), False

def trim_reads(cfg: PipelineConfig, tmp: Path) -> TrimOutputs:
    trim_galore = which_or_die("trim_galore")
    out = tmp / "trim_galore"; out.mkdir(parents=True, exist_ok=True)
    INFO("Trimming reads with Trim Galore (paired + retain unpaired)")
    cmd = [
        trim_galore, "--paired", "--retain_unpaired",
        "--length", str(cfg.length), "--phred33",
        "-o", str(out),
    ]
    if cfg.fastqc:
        cmd.append("--fastqc")
    if cfg.threads and cfg.threads > 1:
        cmd += ["--cores", str(cfg.threads)]
    # add inputs ONCE, then run
    cmd += [cfg.r1, cfg.r2]
    run_cmd(cmd)

    def pick_one(patterns):
        for pat in patterns:
            hits = sorted(out.glob(pat))
            if hits:
                return hits[0]
        return None

    val1 = pick_one(["*_val_1.fq.gz", "*_val_1.fastq.gz"])
    val2 = pick_one(["*_val_2.fq.gz", "*_val_2.fastq.gz"])
    unp1 = pick_one(["*_unpaired_1.fq.gz", "*_unpaired_1.fastq.gz"])
    unp2 = pick_one(["*_unpaired_2.fq.gz", "*_unpaired_2.fastq.gz"])

    if not val1 or not val2:
        listing = "\n".join(sorted(p.name for p in out.iterdir()))
        raise SystemExit(
            f"ERROR: Trim Galore finished but expected outputs not found in {out}.\n"
            f"Looked for *_val_1/2.(fastq|fq).gz\n"
            f"Files present:\n{listing}"
        )
    return TrimOutputs(out, val1, val2, unp1, unp2)

def map_paired(cfg: PipelineConfig, val1: Path, val2: Path, tmp: Path) -> Path:
    bwa      = which_or_die("bwa")
    samtools = which_or_die("samtools")
    INFO("Mapping paired reads")
    pair_bam = tmp / "paired.bam"
    rg = f"@RG\\tID:{cfg.sample}\\tSM:{cfg.sample}\\tPL:illumina"
    cmd = (
        f"{bwa} mem -t {cfg.threads} -K 10000000 -c 100 -R '{rg}' -M -T 50 "
        f"{shlex.quote(cfg.ref)} {shlex.quote(str(val1))} {shlex.quote(str(val2))} | "
        f"{samtools} sort -@ {cfg.threads} -o {shlex.quote(str(pair_bam))} -"
    )
    run_cmd(cmd, shell=True)
    return pair_bam

def map_singles(cfg: PipelineConfig, unp1: Optional[Path], unp2: Optional[Path], tmp: Path) -> Optional[Path]:
    if cfg.skip_singles:
        INFO("Skipping singles mapping as requested")
        return None
    bwa      = which_or_die("bwa")
    samtools = which_or_die("samtools")
    INFO("Preparing unpaired reads (if any)")
    TU = tmp / "singles.fq.gz"
    gz_cat_to(TU, unp1, unp2)
    if not TU.exists() or TU.stat().st_size == 0:
        INFO("No unpaired reads to map (skipping)")
        return None
    INFO("Mapping unpaired reads")
    single_bam = tmp / "single.bam"
    rg = f"@RG\\tID:{cfg.sample}\\tSM:{cfg.sample}\\tPL:illumina"
    cmd = (
        f"{bwa} mem -t {cfg.threads} -K 10000000 -c 100 -R '{rg}' -M -T 50 "
        f"{shlex.quote(cfg.ref)} {shlex.quote(str(TU))} | "
        f"{samtools} sort -@ {cfg.threads} -o {shlex.quote(str(single_bam))} -"
    )
    run_cmd(cmd, shell=True)
    return single_bam

def merge_bams(pair_bam: Path, single_bam: Optional[Path], tmp: Path, threads: int = 1) -> Path:
    samtools = which_or_die("samtools")
    INFO("Merging alignments")
    unsort = tmp / "merged.unsort.bam"
    if single_bam:
        run_cmd([samtools, "merge", "-@", str(threads), "-f", str(unsort), str(pair_bam), str(single_bam)])
    else:
        shutil.copyfile(pair_bam, unsort)
    return unsort

def fix_sort_markdup(cfg: PipelineConfig, unsort_bam: Path) -> str:
    samtools = which_or_die("samtools")
    INFO("Fixmate, sort, mark duplicates")
    out_bam = str(Path(cfg.out_bam or f"{cfg.sample}.bam").resolve())
    chain = (
        f"{samtools} sort -m {cfg.mem} -n -@ {cfg.threads} {shlex.quote(str(unsort_bam))} | "
        f"{samtools} fixmate -@ {cfg.threads} -m - - | "
        f"{samtools} sort -m {cfg.mem} -@ {cfg.threads} - | "
        f"{samtools} markdup -@ {cfg.threads} - {shlex.quote(out_bam)}"
    )
    run_cmd(chain, shell=True)
    return out_bam

def index_bam(bam_path: str, threads: int) -> None:
    samtools = which_or_die("samtools")
    INFO("Indexing BAM")
    run_cmd([samtools, "index", "-@", str(threads), bam_path])

def dump_headers_like_logs(bam_path: str, n: int = 2) -> None:
    samtools = which_or_die("samtools")
    for _ in range(n):
        hdr = str(uuid4())
        DEBUG(f"Running command: samtools view -H {bam_path} > {hdr: <82} utils.py:480")
        with open(hdr, "wb") as w:
            sp.run([samtools, "view", "-H", bam_path], check=True, stdout=w)

def make_reference_indexes(ref: str) -> None:
    samtools = which_or_die("samtools")
    gatk     = shutil.which("gatk")  # may be None if not calling
    refp = Path(ref)
    fai = refp.with_suffix(refp.suffix + ".fai")
    dct = refp.with_suffix(".dict")
    if not fai.exists():
        INFO("Creating FASTA index (.fai)")
        run_cmd([samtools, "faidx", str(refp)])
    if gatk and not dct.exists():
        INFO("Creating sequence dictionary (.dict)")
        run_cmd([gatk, "CreateSequenceDictionary", "-R", str(refp), "-O", str(dct)])

def run_haplotypecaller(cfg: PipelineConfig, bam_path: str) -> str:
    gatk = which_or_die("gatk")
    ref  = cfg.ref
    make_reference_indexes(ref)
    out = Path(f"{cfg.sample}.raw.g.vcf.gz" if cfg.emit_gvcf else f"{cfg.sample}.raw.vcf.gz")
    cmd = [
        gatk, "HaplotypeCaller",
        "-R", ref, "-I", bam_path,
        "--native-pair-hmm-threads", str(cfg.threads),
        "--sample-ploidy", str(cfg.sample_ploidy),
        "-O", str(out),
    ]
    if cfg.emit_gvcf:
        cmd += ["-ERC", "GVCF"]
    if cfg.intervals:
        require_files(cfg.intervals)
        cmd += ["-L", cfg.intervals]
    if cfg.gatk_extra:
        cmd += shlex.split(cfg.gatk_extra)
    INFO("Running GATK HaplotypeCaller")
    run_cmd(cmd)
    return str(out)

def normalize_vcf(vcf_path: str, ref: str, threads: int, split_multiallelics: bool = True) -> str:
    bcftools = which_or_die("bcftools")
    p = Path(vcf_path)
    name = p.name
    if name.endswith(".raw.vcf.gz"):
        out_name = name.replace(".raw.vcf.gz", ".norm.vcf.gz")
    elif name.endswith(".vcf.gz"):
        out_name = name[:-7] + ".norm.vcf.gz"
    else:
        out_name = name + ".norm.vcf.gz"
    out = str(p.with_name(out_name))

    parts = [bcftools, "norm", "--threads", str(threads), "-f", ref]
    if split_multiallelics:
        parts += ["-m", "-"]
    parts += [vcf_path, "-Oz", "-o", out]

    INFO("Normalizing VCF with bcftools norm")
    run_cmd(parts)
    run_cmd([bcftools, "index", "-t", out])
    return out

# ---------------- snpEff (for HGVS/cDNA matching) ----------------
def run_snpeff(vcf_in: str, db: str, threads: int) -> str:
    snpeff   = which_or_die("snpEff")
    bcftools = which_or_die("bcftools")
    out = str(Path(vcf_in).with_name(Path(vcf_in).stem + ".ann.vcf.gz"))
    INFO(f"Annotating with snpEff ({db})")
    # Use bcftools to bgzip+index in one go
    cmd = f"{snpeff} -hgvs {shlex.quote(db)} {shlex.quote(vcf_in)} | {bcftools} view -Oz -o {shlex.quote(out)} -"
    run_cmd(cmd, shell=True)
    run_cmd([bcftools, "index", "-t", out])
    return out

# ---------------- VCF masking & DR annotation ----------------
def apply_mask_bed(vcf_in: str, mask_bed: str, threads: int) -> str:
    bcftools = which_or_die("bcftools")
    require_files(vcf_in, mask_bed)
    p = Path(vcf_in)
    out = str(p.with_name(p.name.replace(".vcf.gz", ".masked.vcf.gz")))
    INFO(f"Applying mask BED to VCF: excluding regions in {mask_bed}")
    run_cmd([bcftools, "view", "--threads", str(threads), "-T", f"^{mask_bed}", vcf_in, "-Oz", "-o", out])
    run_cmd([bcftools, "index", "-t", out])
    return out

def _parse_format(fmt: str, sample: str) -> Dict[str,str]:
    keys = fmt.split(":")
    vals = sample.split(":")
    return {k: (vals[i] if i < len(vals) else "") for i,k in enumerate(keys)}

def _float_or_none(x: str) -> Optional[float]:
    try: return float(x)
    except Exception: return None

RV_RE    = re.compile(r"^Rv\d{4}[A-Za-z]?$")
CDNA_RE  = re.compile(r"^c\.[0-9]+[ACGT]>[ACGT]$")

def _load_dr_json(json_path: str):
    """
    Flex loader:
      - genomic entries -> nt_lut[(chrom,pos,ref,alt)] = [records]
      - cDNA entries    -> cdna_lut[(gene, 'c.xxxN>M')] = [records]
    """
    require_files(json_path)
    with open(json_path) as f:
        data = json.load(f)

    def norm_conf(x):
        x = (str(x) if x is not None else "").strip().lower()
        if x in {"high","strong","certain"}: return "high"
        if x in {"moderate","medium"}: return "moderate"
        return "low"

    nt_lut: DefaultDict[Tuple[str,int,str,str], List[Dict[str,Any]]] = defaultdict(list)
    cdna_lut: DefaultDict[Tuple[str,str], List[Dict[str,Any]]]       = defaultdict(list)

    def push_nt(rec):
        try:
            chrom = str(rec.get("chrom","NC_000962.3"))
            pos   = int(rec["pos"])
            ref   = str(rec["ref"])
            alt   = str(rec["alt"])
            nt_lut[(chrom,pos,ref,alt)].append({
                "drug": str(rec.get("drug","")),
                "gene": str(rec.get("gene","")),
                "change": str(rec.get("change","")),
                "confidence": norm_conf(rec.get("confidence","low")),
            })
        except Exception:
            pass

    def push_cdna(gene, cstr, ctx):
        if not (gene and RV_RE.match(gene) and isinstance(cstr,str) and cstr.startswith("c.")):
            return
        cdna_lut[(gene, cstr)].append({
            "drug": str(ctx.get("drug","")),
            "gene": gene,
            "change": cstr,
            "confidence": norm_conf(ctx.get("confidence","low")),
        })

    # recursive walk to harvest both styles
    def walk(x, ctx):
        if isinstance(x, dict):
            ctx = dict(ctx)
            for k,v in x.items():
                lk = k.lower()
                if lk in {"gene","gene_name","locus","locus_tag","id"} and isinstance(v,str) and RV_RE.match(v.strip()):
                    ctx["gene"] = v.strip()
                elif lk in {"drug","drugs"}:
                    ctx["drug"] = v if isinstance(v,str) else ",".join(v) if isinstance(v,list) else str(v)
                elif lk in {"confidence","who_confidence","conf"}:
                    ctx["confidence"] = v

            # genomic style (has pos/ref/alt)
            if {"pos","ref","alt"}.issubset(set(x.keys())):
                push_nt(x)

            # cDNA/HGVS string variants (many catalogs use 'original_mutation' or 'hgvs' etc.)
            for key in ("original_mutation","hgvs","hgvs_c","nucleotide","nucleotide_change","mutation","change"):
                if key in x and isinstance(x[key], str) and x[key].startswith("c."):
                    push_cdna(ctx.get("gene",""), x[key].strip(), ctx)

            for v in x.values(): walk(v, ctx)
        elif isinstance(x, list):
            for v in x: walk(v, ctx)

    walk(data, {})
    INFO(f"Loaded DR JSON: {sum(len(v) for v in nt_lut.values())} genomic entries; "
         f"{sum(len(v) for v in cdna_lut.values())} cDNA entries")
    return nt_lut, cdna_lut

def _read_vcf_records(vcf_path: str):
    require_files(vcf_path)
    if vcf_path.endswith(".vcf.gz"):
        p = sp.Popen(["bcftools","view","-H", vcf_path], stdout=sp.PIPE, text=True)
        try:
            for line in p.stdout:
                yield line.rstrip("\n")
        finally:
            p.stdout.close()
    else:
        with open(vcf_path) as f:
            for line in f:
                if line and not line.startswith("#"):
                    yield line.rstrip("\n")

def annotate_dr_from_vcf(vcf_path: str, catalog: Dict[Tuple[str,int,str,str], List[Dict[str,Any]]],
                         min_af: float = 0.0, min_dp: int = 0) -> Tuple[List[Dict[str,Any]], Dict[str,str]]:
    muts: List[Dict[str,Any]] = []
    hits_by_drug: DefaultDict[str, List[Dict[str,Any]]] = defaultdict(list)

    for line in _read_vcf_records(vcf_path):
        fields = line.split("\t")
        if len(fields) < 8: continue
        chrom, pos_s, _id, ref, alts_s, qual, flt, info = fields[:8]
        fmt    = fields[8] if len(fields) > 8 else ""
        sample = fields[9] if len(fields) > 9 else ""
        pos = int(pos_s)
        alt_list = alts_s.split(",")

        fmt_map = _parse_format(fmt, sample) if sample else {}
        dp = int(fmt_map.get("DP","0") or 0)
        ad = fmt_map.get("AD","")
        ref_count, alt_counts = None, []
        if ad:
            parts = [p for p in ad.split(",") if p != ""]
            if parts:
                ref_count = _float_or_none(parts[0])
                alt_counts = [ _float_or_none(x) for x in parts[1:1+len(alt_list)] ]

        for i, alt in enumerate(alt_list):
            key = (chrom, pos, ref, alt)
            if key not in catalog:
                continue
            ac = alt_counts[i] if i < len(alt_counts) else None
            af = (ac / ((ref_count or 0.0) + ac)) if (ac is not None and (ref_count or 0.0)+ac>0) else None
            if af is None and "AF=" in info:
                try:
                    af = float(info.split("AF=")[1].split(";")[0].split(",")[i])
                except Exception:
                    pass
            if min_dp and dp < min_dp: 
                continue
            if min_af and (af is None or af < min_af):
                continue

            for rec in catalog[key]:
                m = {
                    "drug": rec["drug"],
                    "gene": rec.get("gene",""),
                    "chrom": chrom, "pos": pos, "ref": ref, "alt": alt,
                    "change": rec.get("change",""),
                    "confidence": rec["confidence"],
                    "DP": dp,
                    "AF": round(af,4) if af is not None else None,
                    "filter": flt if flt not in (".","PASS") else "PASS"
                }
                muts.append(m)
                hits_by_drug[rec["drug"]].append(m)

    per_drug: Dict[str,str] = {}
    for drug, hits in hits_by_drug.items():
        confs = {h["confidence"] for h in hits}
        if "high" in confs or "moderate" in confs: per_drug[drug] = "R"
        elif "low" in confs: per_drug[drug] = "U"
        else: per_drug[drug] = "U"
    return muts, per_drug

def annotate_dr_from_ann_vcf(vcf_path: str,
                             catalog_cdna: Dict[Tuple[str,str], List[Dict[str,Any]]],
                             min_af: float = 0.0, min_dp: int = 0):
    muts: List[Dict[str,Any]] = []
    hits_by_drug: DefaultDict[str, List[Dict[str,Any]]] = defaultdict(list)

    for line in _read_vcf_records(vcf_path):
        f = line.split("\t")
        if len(f) < 8: continue
        chrom, pos_s, _id, ref, alts_s, qual, flt, info = f[:8]
        fmt    = f[8] if len(f) > 8 else ""
        sample = f[9] if len(f) > 9 else ""
        pos = int(pos_s)
        alt_list = alts_s.split(",")

        fmt_map = _parse_format(fmt, sample) if sample else {}
        dp = int(fmt_map.get("DP","0") or 0)
        ad = fmt_map.get("AD","")
        ref_count, alt_counts = None, []
        if ad:
            parts = [p for p in ad.split(",") if p != ""]
            if parts:
                ref_count = _float_or_none(parts[0])
                alt_counts = [_float_or_none(x) for x in parts[1:1+len(alt_list)]]

        # pull ANN entries
        ann_items = []
        for kv in info.split(";"):
            if kv.startswith("ANN="):
                ann_items = kv[4:].split(",")
                break

        # ANN fields: Allele|Annotation|Impact|Gene_Name|Gene_ID|Feature_Type|Feature_ID|BioType|Rank|HGVS.c|HGVS.p|...
        for ann in ann_items:
            cols = ann.split("|")
            if len(cols) < 11: 
                continue
            gene_name = cols[3] or ""
            gene_id   = cols[4] or ""
            cdna      = cols[9] or ""
            if not cdna.startswith("c."):
                continue
            gene_key = gene_id if RV_RE.match(gene_id or "") else (gene_name if RV_RE.match(gene_name or "") else None)
            if not gene_key:
                continue
            key = (gene_key, cdna)
            if key not in catalog_cdna:
                continue

            # thresholds
            if min_dp and dp < min_dp:
                continue

            af = None
            if alt_counts:
                acs = [a for a in alt_counts if a is not None]
                if acs and (ref_count or 0) + max(acs) > 0:
                    af = max(acs) / ((ref_count or 0) + max(acs))
            if af is None and "AF=" in info:
                try:
                    afs = [float(x) for x in info.split("AF=")[1].split(";")[0].split(",")]
                    af = max(afs) if afs else None
                except Exception:
                    pass
            if min_af and (af is None or af < min_af):
                continue

            for rec in catalog_cdna[key]:
                m = {
                    "drug": rec["drug"],
                    "gene": gene_key,
                    "chrom": chrom, "pos": pos, "ref": ref, "alt": alts_s,
                    "change": cdna,
                    "confidence": rec["confidence"],
                    "DP": dp,
                    "AF": round(af,4) if af is not None else None,
                    "filter": flt if flt not in (".","PASS") else "PASS",
                }
                muts.append(m)
                hits_by_drug[rec["drug"]].append(m)

    per_drug: Dict[str,str] = {}
    for drug, hits in hits_by_drug.items():
        confs = {h["confidence"] for h in hits}
        if "high" in confs or "moderate" in confs: per_drug[drug] = "R"
        elif "low" in confs: per_drug[drug] = "U"
        else: per_drug[drug] = "U"
    return muts, per_drug

def write_dr_outputs(sample: str, outdir: str, mutations: List[Dict[str,Any]], per_drug: Dict[str,str]) -> Dict[str,str]:
    outp = {}
    outdir_p = Path(outdir)
    outdir_p.mkdir(parents=True, exist_ok=True)

    muts_tsv = outdir_p / f"{sample}.dr_mutations.tsv"
    with open(muts_tsv, "w", newline="") as w:
        cols = ["drug","gene","chrom","pos","ref","alt","change","confidence","DP","AF","filter"]
        w.write("\t".join(cols) + "\n")
        for m in sorted(mutations, key=lambda x: (x["drug"], x["gene"], x["chrom"], x["pos"])):
            w.write("\t".join(str(m.get(c,"")) for c in cols) + "\n")
    outp["dr_mutations"] = str(muts_tsv)

    summ_tsv = outdir_p / f"{sample}.dr_summary.tsv"
    all_drugs = sorted(per_drug.keys())
    with open(summ_tsv, "w", newline="") as w:
        w.write("sample\tdrug\tprediction\tn_mutations\n")
        for d in all_drugs:
            n = sum(1 for m in mutations if m["drug"] == d)
            w.write(f"{sample}\t{d}\t{per_drug[d]}\t{n}\n")
    outp["dr_summary"] = str(summ_tsv)

    js = {
        "sample": sample,
        "timestamp": time.strftime("%Y-%m-%d %H:%M:%S"),
        "drugs": [
            {"drug": d, "prediction": per_drug[d],
             "mutations": [m for m in mutations if m["drug"] == d]}
            for d in all_drugs
        ],
        "n_total_mutations": len(mutations)
    }
    js_path = outdir_p / f"{sample}.tbprofiler_like.json"
    with open(js_path, "w") as w:
        json.dump(js, w, indent=2)
    outp["dr_json"] = str(js_path)

    INFO(f"DR outputs written: {js_path.name}, {summ_tsv.name}, {muts_tsv.name}")
    return outp

# ---------------- orchestration ----------------
def run_pipeline(cfg: PipelineConfig) -> dict[str, str]:
    require_files(cfg.r1, cfg.r2, cfg.ref)
    tmp, keep_tmp = prepare_tempdir(cfg)
    INFO(f"Reference: {cfg.ref}")
    INFO(f"Temp dir: {tmp}")

    with step("Trimming reads", "fastq.py:39"):
        trims = trim_reads(cfg, tmp)

    with step("Mapping paired reads", "fastq.py:53"):
        pair_bam   = map_paired(cfg, trims.val1, trims.val2, tmp)

    with step("Mapping unpaired reads (if any)", "fastq.py:53"):
        single_bam = map_singles(cfg, trims.unp1, trims.unp2, tmp)

    with step("Merging alignments", "utils.py:480"):
        unsort_bam = merge_bams(pair_bam, single_bam, tmp, cfg.threads)

    with step("Fixmate, sort, mark duplicates", "bam.py:36"):
        final_bam  = fix_sort_markdup(cfg, unsort_bam)

    with step("Indexing BAM", "utils.py:480"):
        index_bam(final_bam, cfg.threads)

    with step("Dumping BAM headers", "utils.py:480"):
        dump_headers_like_logs(final_bam, n=2)

    outputs: Dict[str, str] = {"bam": final_bam, "bai": final_bam + ".bai"}

    vcf_for_annotation: Optional[str] = None
    if cfg.call_variants:
        with step("Running GATK HaplotypeCaller", "bam.py:162"):
            raw_vcf = run_haplotypecaller(cfg, final_bam)
            outputs["vcf_raw"] = raw_vcf
        vcf_for_annotation = raw_vcf

        if cfg.normalize:
            with step("Normalizing VCF with bcftools", "vcf.py:119"):
                norm = normalize_vcf(raw_vcf, cfg.ref, cfg.threads, split_multiallelics=not cfg.emit_gvcf)
                outputs["vcf_norm"] = norm
                vcf_for_annotation = norm

    # Apply mask before annotation (TB-Profiler style)
        if cfg.mask_bed and vcf_for_annotation:
            with step("Masking VCF with BED", "vcf.py:141"):
                vcf_for_annotation = apply_mask_bed(vcf_for_annotation, cfg.mask_bed, cfg.threads)
                outputs["vcf_masked"] = vcf_for_annotation

    # Annotate drug resistance
    if cfg.dr_catalog and vcf_for_annotation:
        with step("Drug resistance annotation", "dr.py:42"):
            nt_lut, cdna_lut = _load_dr_json(cfg.dr_catalog)

            vcf_ann = vcf_for_annotation
            if cdna_lut:
                # add snpEff annotations so we can read HGVS.c & gene
                with step("snpEff annotation", "ann.py:10"):
                    vcf_ann = run_snpeff(vcf_for_annotation, cfg.snpeff_db or "Mycobacterium_tuberculosis_h37rv", cfg.threads)
                    outputs["vcf_ann"] = vcf_ann

            muts_all: List[Dict[str,Any]] = []
            per_drug_all: Dict[str,str] = {}

            if nt_lut:
                m1, d1 = annotate_dr_from_vcf(vcf_ann, nt_lut, min_af=cfg.min_alt_af, min_dp=cfg.min_dp)
                muts_all += m1
                per_drug_all.update(d1)

            if cdna_lut:
                m2, d2 = annotate_dr_from_ann_vcf(vcf_ann, cdna_lut, min_af=cfg.min_alt_af, min_dp=cfg.min_dp)
                muts_all += m2
                for k,v in d2.items():
                    per_drug_all[k] = "R" if (v == "R" or per_drug_all.get(k) == "R") else (per_drug_all.get(k) or v or "U")

            if not muts_all and not per_drug_all:
                WARNING("No DR matches found. If your catalog is cDNA-based, ensure snpEff DB is installed.")

            dr_paths = write_dr_outputs(cfg.sample, Path(outputs["bam"]).parent, muts_all, per_drug_all)
            outputs.update(dr_paths)

    if keep_tmp:
        INFO(f"Temp dir kept: {tmp}")
    else:
        INFO("Temp dir was ephemeral and has been removed")

    INFO(f"Done. Final BAM: {final_bam}")
    return outputs

# ---------------- CLI ----------------
def _build_argparser():
    import argparse
    p = argparse.ArgumentParser(description="Trim->Map->HC->norm + optional DR annotation (TB-Profiler-style).")
    # IO
    p.add_argument("--r1", required=True, help="R1 FASTQ(.gz)")
    p.add_argument("--r2", required=True, help="R2 FASTQ(.gz)")
    p.add_argument("--ref", required=True, help="Reference FASTA")
    p.add_argument("--sample", required=True, help="Sample name")
    p.add_argument("--out-bam", default=None, help="Output BAM path (optional)")
    p.add_argument("--tmpdir", default=None, help="Keep intermediates under this directory")
    # mapping/trimming
    p.add_argument("--length", type=int, default=36, help="Minimum read length after trimming")
    p.add_argument("--threads", type=int, default=1)
    p.add_argument("--mem", default="768M")
    p.add_argument("--fastqc", action="store_true")
    p.add_argument("--skip-singles", action="store_true")
    # calling
    p.add_argument("--call-variants", action="store_true")
    p.add_argument("--emit-gvcf", action="store_true")
    p.add_argument("--intervals", default=None)
    p.add_argument("--normalize", action="store_true")
    p.add_argument("--gatk-extra", default=None)
    # DR
    p.add_argument("--dr-catalog", default=None, help="Drug-resistance JSON (e.g., tbdb.dr.json)")
    p.add_argument("--mask-bed", default=None, help="BED to exclude regions before DR annotation")
    p.add_argument("--min-alt-af", type=float, default=0.0, help="Min ALT AF to count")
    p.add_argument("--min-dp", type=int, default=0, help="Min depth to count")
    p.add_argument("--snpeff-db", default="Mycobacterium_tuberculosis_h37rv",
                   help="snpEff DB name (default: Mycobacterium_tuberculosis_h37rv)")
    # Logging
    p.add_argument("--log-level", default="INFO", choices=list(_LEVELS.keys()))
    p.add_argument("--log-file", default=None)
    return p

def main():
    args = _build_argparser().parse_args()
    set_log_level(args.log_level)
    set_log_file(args.log_file)
    try:
        cfg = PipelineConfig(
            r1=args.r1, r2=args.r2, ref=args.ref, sample=args.sample,
            out_bam=args.out_bam, tmpdir=args.tmpdir, threads=args.threads,
            mem=args.mem, fastqc=args.fastqc, skip_singles=args.skip_singles,
            call_variants=args.call_variants, emit_gvcf=args.emit_gvcf,
            intervals=args.intervals, normalize=args.normalize, gatk_extra=args.gatk_extra,
            dr_catalog=args.dr_catalog, min_alt_af=args.min_alt_af, min_dp=args.min_dp,
            mask_bed=args.mask_bed, snpeff_db=args.snpeff_db
        )
        outs = run_pipeline(cfg)
        INFO("Outputs: " + " ".join(f"{k}={v}" for k,v in outs.items()))
    finally:
        close_logs()

if __name__ == "__main__":
    main()

