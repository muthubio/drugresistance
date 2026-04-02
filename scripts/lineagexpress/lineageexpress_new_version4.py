#!/usr/bin/env python3
import os, subprocess, argparse, logging, time, re
from collections import defaultdict
from pathlib import Path
import pandas as pd


# --------------------------- Logging ---------------------------
logging.basicConfig(format='%(asctime)s - %(levelname)s - %(message)s',
                    level=logging.INFO)


# --------------------------- Defaults --------------------------
BASE_DIR = Path(__file__).resolve().parent
DATA_DIR = BASE_DIR / "data"
DEFAULT_REF       = str(DATA_DIR / "h37rv.fa")
DEFAULT_SNP_FILE  = str(DATA_DIR / "lineage_snp_updated_au13.tsv")
DEFAULT_BED_FILE  = str(DATA_DIR / "targeted_modified_regions_au13.bed")
DEFAULT_SNPEFF_DB = "Mycobacterium_tuberculosis_h37rv"


SAMPLE_INPUT_DIR = BASE_DIR / "sample_input"
SAMPLE_DATA_DIR  = BASE_DIR / "sample_data"


# ------------------------ Shell runner -------------------------
def run_cmd(cmd: str) -> bool:
    logging.info(f"Running: {cmd}")
    try:
        subprocess.run(cmd, check=True, shell=True)
        return True
    except subprocess.CalledProcessError as e:
        logging.error(f"Command failed: {cmd}\n{e}")
        return False


# -------------------- Flexible input resolver ------------------
def _find_in_dir(dirpath: Path, sid: str):
    """Search common paired/single patterns inside a directory."""
    pats = [
        (f"{sid}_1.fastq.gz",  f"{sid}_2.fastq.gz"),
        (f"{sid}_R1.fastq.gz", f"{sid}_R2.fastq.gz"),
        (f"{sid}.1.fastq.gz",  f"{sid}.2.fastq.gz"),
        (f"{sid}_1.fastq",     f"{sid}_2.fastq"),
        (f"{sid}_R1.fastq",    f"{sid}_R2.fastq"),
        (f"{sid}.1.fastq",     f"{sid}.2.fastq"),
    ]
    for r1, r2 in pats:
        r1p, r2p = dirpath / r1, dirpath / r2
        if r1p.exists() and r2p.exists():
            return "paired", [str(r1p), str(r2p)]
    for ext in (".fastq.gz",".fq.gz",".fastq",".fq"):
        se = dirpath / f"{sid}{ext}"
        if se.exists():
            return "single", [str(se)]
    try:
        fqs = sorted([p for p in dirpath.iterdir()
                      if p.is_file() and (p.suffix in {".fastq",".fq"} or p.name.endswith(".fastq.gz") or p.name.endswith(".fq.gz"))])
    except FileNotFoundError:
        fqs = []
    if len(fqs) == 2:
        return "paired", [str(fqs[0]), str(fqs[1])]
    if len(fqs) == 1:
        return "single", [str(fqs[0])]
    return "unknown", []


def resolve_inputs(sample_spec: str):
    """
    Accepts:
      - 'R1.fastq.gz,R2.fastq.gz'
      - file path (fastq/fq/bam/vcf)
      - directory path (scan for reads)
      - bare sample ID (search CWD, sample_input/, sample_data/)
    Returns: (kind, [paths]) with kind in {'paired','single','bam','vcf','unknown'}
    """
    s = sample_spec.strip()
    if "," in s:
        r1s, r2s = [x.strip() for x in s.split(",", 1)]
        r1p, r2p = Path(r1s).expanduser(), Path(r2s).expanduser()
        if r1p.exists() and r2p.exists():
            return "paired", [str(r1p), str(r2p)]
        return "unknown", []


    p = Path(s).expanduser()


    if p.is_file():
        low = p.name.lower()
        if low.endswith((".vcf", ".vcf.gz")): return "vcf", [str(p)]
        if low.endswith(".bam"):               return "bam", [str(p)]
        if low.endswith((".fastq",".fq",".fastq.gz",".fq.gz")):
            return "single", [str(p)]
        return "unknown", []


    if p.is_dir():
        kind, files = _find_in_dir(p, p.name)
        if kind != "unknown": return kind, files
        return _find_in_dir(p, p.name)


    sid = s
    for base in [Path.cwd(), SAMPLE_INPUT_DIR, SAMPLE_DATA_DIR]:
        kind, files = _find_in_dir(base, sid)
        if kind != "unknown":
            return kind, files


    return "unknown", []


# ---------------------- Lineage utilities ----------------------
def load_lineage_snp_file(snp_file):
    lineage_snps = defaultdict(set)
    with open(snp_file, 'r') as f:
        for line in f:
            if line.strip() and not line.startswith('#'):
                parts = line.strip().split("\t")
                if len(parts) == 3:
                    lineage, rvlocus, pos = parts
                elif len(parts) == 2:
                    lineage, pos = parts
                    rvlocus = "NC_000962.3"
                else:
                    logging.warning(f"Skipping malformed line: {line.strip()}")
                    continue
                try:
                    lineage_snps[lineage].add((rvlocus, int(pos)))
                except ValueError:
                    logging.warning(f"Bad position in line: {line.strip()}")
    logging.info(f"Loaded SNP reference data for {len(lineage_snps)} lineages from {snp_file}")
    return lineage_snps


def pct(matched, total): return round((matched / total) * 100, 2) if total > 0 else 0.0


# -------- Contig alias bridge: NC_000962.3 <-> Chromosome -----
CONTIG_ALIASES = {"Chromosome", "NC_000962.3"}


def contig_alias_to_nc(chrom: str) -> str:
    """Normalize to 'NC_000962.3' for lineage comparison."""
    return "NC_000962.3" if chrom in CONTIG_ALIASES else chrom


def rewrite_bed_contig(bed_in: str, bed_out: str, target: str = "NC_000962.3"):
    """Create BED copy whose contig matches `target` if an alias is seen."""
    with open(bed_in) as fi, open(bed_out, "w") as fo:
        for ln in fi:
            if not ln.strip() or ln.startswith("#"):
                fo.write(ln); continue
            parts = ln.rstrip("\n").split("\t")
            if parts:
                parts[0] = target if parts[0] in CONTIG_ALIASES else parts[0]
            fo.write("\t".join(parts) + "\n")


def lineage_snps_to_target(lineage_snps: dict, target: str = "NC_000962.3"):
    converted = {}
    for lin, sset in lineage_snps.items():
        converted[lin] = {(target if c in CONTIG_ALIASES else c, pos) for (c,pos) in sset}
    return converted


def write_chrmap(tmp_map_path: str, to_chr: bool = True):
    """bcftools rename map."""
    with open(tmp_map_path, "w") as fh:
        if to_chr:
            fh.write("NC_000962.3\tChromosome\n")
        else:
            fh.write("Chromosome\tNC_000962.3\n")


# --------------- snpEff parsing & matching utils ---------------
AA3_TO_1 = {
    "Ala":"A","Arg":"R","Asn":"N","Asp":"D","Cys":"C","Gln":"Q","Glu":"E","Gly":"G","His":"H",
    "Ile":"I","Leu":"L","Lys":"K","Met":"M","Phe":"F","Pro":"P","Ser":"S","Thr":"T","Trp":"W",
    "Tyr":"Y","Val":"V","Ter":"*","Stop":"*"
}
def _n(s): return ("" if pd.isna(s) else str(s)).strip()
def _g(s): return _n(s).lower()
def _rv(s): return _n(s).lower()
def _c(s): return _n(s).replace(" ", "")


def to_short_protein(hgvs_p: str) -> str:
    s = _n(hgvs_p)
    if not s: return ""
    if s.startswith("p."):
        m = re.fullmatch(r"p\.([A-Z][a-z]{2})(\d+)([A-Z][a-z]{2}|\*)", s)
        if m:
            a1,pos,a2 = m.groups()
            a1 = AA3_TO_1.get(a1, "")
            a2 = "*" if a2 in ("*", "Ter", "Stop") else AA3_TO_1.get(a2, "")
            if a1 and a2: return f"{a1}{pos}{a2}"
    m = re.search(r"([A-Z\*]\d+[A-Z\*])", s)
    return m.group(1) if m else ""


def parse_snpeff_vcf_to_df(ann_vcf_path: str) -> pd.DataFrame:
    rows = []
    with open(ann_vcf_path, "r", encoding="utf-8", errors="ignore") as fh:
        for line in fh:
            if not line or line.startswith("#"): continue
            parts = line.rstrip("\n").split("\t", 8)
            if len(parts) < 8: continue
            chrom, pos, _id, ref, alt, _qual, _filter, info = parts[:8]
            ann_val = None
            for kv in info.split(";"):
                if kv.startswith("ANN="):
                    ann_val = kv[4:]; break
            if not ann_val: continue
            first_ann = ann_val.split(",", 1)[0]
            cols = first_ann.split("|")
            if len(cols) < 16: cols += [""] * (16 - len(cols))
            rows.append([
                int(pos), ref, alt,
                cols[1], cols[2], cols[3], cols[4], cols[9], cols[10],
            ])
    return pd.DataFrame(rows, columns=[
        "Position","REF","ALT","Annotation","Impact","Gene","Rv_locus","NA_change","AA_change"
    ])


# ---- Confidence predicate: permissive and robust for WHO CSV ---
def conf_is_resistant(conf: str) -> bool:
    """
    True for WHO-style positives; False for 'Not assoc w R' etc.
    """
    if not isinstance(conf, str):
        return False
    s = conf.strip().lower()
    if "not assoc" in s or "not associated" in s:
        return False
    return ("assoc" in s) or ("associated" in s) or ("confers" in s) or ("causes" in s) or ("likely" in s)


# ------------------- Catalogue matching ------------------------
def match_catalogue(ann_csv: str, mutations_csv: str,
                    out_all: str, out_res: str):
    # --- annotated table ---
    ann = pd.read_csv(ann_csv)
    ann["Gene_norm"] = ann.get("Gene", "").map(_g)
    ann["Rv_norm"]   = ann.get("Rv_locus", pd.Series([""]*len(ann))).map(_rv)
    ann["cDNA_norm"] = ann.get("NA_change", pd.Series([""]*len(ann))).map(_c)
    ann["AA_short"]  = ann.get("AA_change", pd.Series([""]*len(ann))).map(to_short_protein)
    ann = ann[(ann["Gene_norm"] != "") | (ann["Rv_norm"] != "")].copy()


    # --- catalogue ---
    muts = pd.read_csv(mutations_csv)
    if "Gene" not in muts.columns: muts["Gene"] = ""
    if "Mutation" not in muts.columns: muts["Mutation"] = ""
    if "original_mutation" not in muts.columns: muts["original_mutation"] = ""
    muts["Gene_norm"] = muts["Gene"].map(_g)
    muts["mut_pick"]  = muts["Mutation"].fillna(muts["original_mutation"])
    muts["cDNA_norm"] = muts["mut_pick"].map(_c)
    muts["AA_short"]  = muts["mut_pick"].map(to_short_protein)
    muts["Gene_or_rv_norm"] = muts["Gene"].map(_rv)


    matches = []


    # T1: Gene + cDNA
    if not ann.empty and not muts.empty:
        t1 = ann.merge(muts, how="inner",
                       left_on=["Gene_norm","cDNA_norm"],
                       right_on=["Gene_norm","cDNA_norm"],
                       suffixes=("_ann","_cat"))
        if not t1.empty: t1 = t1.copy(); t1["_tier"] = "T1_gene+cDNA"; matches.append(t1)


    # T2: Rv + cDNA
    if not ann.empty and ann["Rv_norm"].ne("").any():
        t2 = ann.merge(muts, how="inner",
                       left_on=["Rv_norm","cDNA_norm"],
                       right_on=["Gene_or_rv_norm","cDNA_norm"],
                       suffixes=("_ann","_cat"))
        if not t2.empty: t2 = t2.copy(); t2["_tier"] = "T2_rv+cDNA"; matches.append(t2)


    # T3: Gene/Rv + protein
    ann_p = ann[ann["AA_short"]!=""].copy()
    muts_p = muts[muts["AA_short"]!=""].copy()
    if not ann_p.empty and not muts_p.empty:
        t3 = ann_p.merge(muts_p, how="inner",
                         left_on=["Gene_norm","AA_short"],
                         right_on=["Gene_norm","AA_short"],
                         suffixes=("_ann","_cat"))
        if not t3.empty: t3 = t3.copy(); t3["_tier"] = "T3_gene+AA"; matches.append(t3)
        if ann["Rv_norm"].ne("").any():
            t3b = ann_p.merge(muts_p, how="inner",
                              left_on=["Rv_norm","AA_short"],
                              right_on=["Gene_or_rv_norm","AA_short"],
                              suffixes=("_ann","_cat"))
            if not t3b.empty: t3b = t3b.copy(); t3b["_tier"] = "T3b_rv+AA"; matches.append(t3b)


    if not matches:
        pd.DataFrame().to_csv(out_all, index=False)
        pd.DataFrame().to_csv(out_res, index=False)
        logging.info("No catalogue matches found.")
        return pd.DataFrame(), pd.DataFrame()


    merged = pd.concat(matches, ignore_index=True, sort=False)


    # Ensure a plain 'Gene' column exists
    if "Gene" not in merged.columns:
        if "Gene_ann" in merged.columns: merged["Gene"] = merged["Gene_ann"]
        elif "Gene_cat" in merged.columns: merged["Gene"] = merged["Gene_cat"]
        else: merged["Gene"] = ""


    keep_context = [c for c in [
        "Position","REF","ALT","Annotation","Impact","Gene","Rv_locus",
        "NA_change","AA_change","_tier"
    ] if c in merged.columns]
    keep_catalog = [c for c in [
        "Gene","Mutation","original_mutation","type","drug","confidence","source","comment"
    ] if c in merged.columns]


    cols_all = list(dict.fromkeys(keep_context + keep_catalog))
    out_all_df = merged[cols_all].copy()


    mutkey = out_all_df.get("Mutation").fillna(out_all_df.get("original_mutation"))
    dedup_key = pd.DataFrame({
        "Gene": out_all_df.get("Gene"),
        "MutKey": mutkey,
        "drug": out_all_df.get("drug"),
        "tier": out_all_df.get("_tier"),
        "Position": out_all_df.get("Position"),
    })
    out_all_df = out_all_df.loc[~dedup_key.duplicated()].reset_index(drop=True)


    if "Position" in out_all_df.columns:
        try: out_all_df["Position"] = pd.to_numeric(out_all_df["Position"], errors="coerce")
        except Exception: pass


    sort_order = [c for c in ["Gene","Position","drug","_tier"] if c in out_all_df.columns]
    if sort_order: out_all_df = out_all_df.sort_values(sort_order, kind="mergesort")


    out_all_df.to_csv(out_all, index=False)


    # General results: all matches (not filtered by resistance)
    out_res_df = out_all_df.copy()


    logging.info(f"Total matches found: {len(out_res_df)}")


    out_res_df.to_csv(out_res, index=False)


    logging.info(f"wrote {out_all}")
    logging.info(f"wrote {out_res}")
    return out_all_df, out_res_df


# ---------------- Drug summary (per-sample) --------------------
def write_drug_summary(matched_res_path: str, out_csv: str):
    try:
        df = pd.read_csv(matched_res_path)
    except Exception:
        logging.info("No matched_mutations.csv to summarize.")
        return
    if df.empty:
        logging.info("matched_mutations.csv is empty; creating sensitive summary.")
        sensitive_summary = pd.DataFrame({
            "drug": ["pan-susceptible"],
            "Gene": ["N/A"],
            "Mutation": ["No resistance mutations"],
            "confidence": ["Sensitive"],
            "n": [1]
        })
        sensitive_summary.to_csv(out_csv, index=False)
        logging.info(f"wrote sensitive {out_csv}")
        return
    for c in ["drug","Gene","Mutation","original_mutation","confidence"]:
        if c not in df.columns: df[c] = ""
    df["MutKey"] = df["Mutation"].fillna(df["original_mutation"])
    summary = (df.groupby(["drug","Gene","MutKey","confidence"], dropna=False)
               .size().reset_index(name="n")
               .sort_values(["drug","Gene","MutKey"]))
    summary = summary.rename(columns={"MutKey":"Mutation"})
    summary.to_csv(out_csv, index=False)
    logging.info(f"wrote {out_csv}")


# --------------- WHO phenotype prediction ---------------------
def _norm(s): return ("" if pd.isna(s) else str(s)).strip().lower()
RIF_NAMES = {"rifampicin","rifampin"}
INH_NAMES = {"isoniazid","inh"}
FQ_NAMES  = {"levofloxacin","moxifloxacin","ofloxacin","gatifloxacin"}
GA_EXTRA  = {"bedaquiline","linezolid"}


def _present(df, name_set): 
    if "confidence" not in df.columns:
        return False
    resistant_df = df[df["confidence"].map(conf_is_resistant)]
    return resistant_df["drug"].map(_norm).isin(name_set).any()

def _evidence_rows(df, name_set, cols=("drug","Gene","Mutation","original_mutation","confidence")):
    if "confidence" not in df.columns:
        return pd.DataFrame()
    resistant_df = df[df["confidence"].map(conf_is_resistant)].copy()
    cols = [c for c in cols if c in resistant_df.columns]
    return (resistant_df[resistant_df["drug"].map(_norm).isin(name_set)][cols]
            .fillna("").drop_duplicates().reset_index(drop=True))


def predict_final_category(matched_res_csv: str, out_txt: str, out_csv: str) -> str:
    if not Path(matched_res_csv).exists():
        logging.warning(f"{matched_res_csv} not found")
        return "No data"
    df = pd.read_csv(matched_res_csv)
    if df.empty:
        cat = "Sensitive Strain"
        Path(out_txt).write_text(cat + "\n", encoding="utf-8")
        pd.DataFrame({"category":[cat]}).to_csv(out_csv, index=False)
        logging.info(cat); return cat
    if "drug" not in df.columns:
        raise SystemExit("matched_mutations.csv has no 'drug' column")


    has_rif = _present(df, RIF_NAMES)
    has_inh = _present(df, INH_NAMES)
    has_fq  = _present(df, FQ_NAMES)
    has_ga  = _present(df, GA_EXTRA)


    if (has_rif or has_inh) and has_fq and has_ga:
        category = "XDR-TB"
    elif (has_rif and has_inh) and has_fq:
        category = "pre-XDR-TB"
    elif has_rif and has_fq:
        category = "pre-XDR-TB"
    elif has_rif and has_inh:
        category = "MDR-TB"
    elif has_rif:
        category = "RR-TB"
    elif has_inh and not has_rif:
        category = "Hr-TB"
    else:
        category = "Sensitive Strain"


    ev_rif = _evidence_rows(df, RIF_NAMES)
    ev_inh = _evidence_rows(df, INH_NAMES)
    ev_fq  = _evidence_rows(df, FQ_NAMES)
    ev_ga  = _evidence_rows(df, GA_EXTRA)


    lines = [f"Final category: {category}\n"]
    def add_block(title, ev):
        if not ev.empty:
            lines.append(f"{title}:")
            for _, r in ev.iterrows():
                mut = r.get("Mutation") or r.get("original_mutation") or ""
                lines.append(f"  - {r.get('drug','')} | {r.get('Gene','')} | {mut} | {r.get('confidence','')}")
            lines.append("")
    add_block("Rifampicin evidence", ev_rif)
    add_block("Isoniazid evidence", ev_inh)
    add_block("Fluoroquinolone evidence", ev_fq)
    add_block("Group A (Bedaquiline/Linezolid) evidence", ev_ga)


    Path(out_txt).write_text("\n".join(lines).rstrip()+"\n", encoding="utf-8")
    pd.DataFrame([{
        "RR_any_rif": bool(has_rif),
        "Hr_any_inh": bool(has_inh),
        "FQ_any": bool(has_fq),
        "GroupA_additional_any": bool(has_ga),
        "category": category
    }]).to_csv(out_csv, index=False)


    logging.info(f"Final category: {category}")
    return category


# -------------------- per-sample processing --------------------
def process_sample(sample_path, ref_genome, output_dir, lineage_references,
                   bed_file, threads=1,
                   bam_override=None, vcf_override=None,
                   snpeff_cmd="snpEff", snpeff_db=DEFAULT_SNPEFF_DB,
                   mutations_csv=None):
    start_time = time.time()
    output_dir = Path(output_dir); output_dir.mkdir(parents=True, exist_ok=True)
    logs_dir = output_dir / "logs"; logs_dir.mkdir(exist_ok=True)


    sample_id = Path(str(sample_path)).name
    sample_id = sample_id.split(",")[0]
    sample_id = Path(sample_id).stem
    sample_dir = output_dir / sample_id; sample_dir.mkdir(exist_ok=True)


    vcf_lineage = None
    vcf_resist  = None
    sorted_bam  = None


    # Resolve inputs (unless overrides provided)
    if vcf_override and Path(vcf_override).exists():
        # Use the same VCF for both flows if user provided one (rare)
        vcf_lineage = str(Path(vcf_override))
        vcf_resist  = str(Path(vcf_override))
        logging.info(f"[{sample_id}] Using provided VCF for both lineage & resistance: {vcf_override}")
    else:
        if bam_override and Path(bam_override).exists():
            sorted_bam = str(Path(bam_override))
            logging.info(f"[{sample_id}] Using provided BAM: {sorted_bam}")
        else:
            kind, inputs = resolve_inputs(str(sample_path))
            if kind == "unknown":
                logging.error(f"[{sample_id}] Could not resolve inputs for '{sample_path}'.")
                return
            if kind == "vcf":
                vcf_lineage = inputs[0]
                vcf_resist  = inputs[0]
                logging.info(f"[{sample_id}] Detected VCF: {vcf_lineage}")
            elif kind == "bam":
                sorted_bam = inputs[0]
                logging.info(f"[{sample_id}] Detected BAM: {sorted_bam}")
            else:
                # FASTQ -> align
                sam_file = sample_dir / f"{sample_id}.sam"
                bwa_log  = logs_dir / f"{sample_id}.bwa.log"
                fastq_str = " ".join(inputs)
                bwa_cmd = (
                    f"bwa mem -M -t {threads} -R '@RG\\tID:{sample_id}\\tSM:{sample_id}\\tPL:ILLUMINA' "
                    f"{ref_genome} {fastq_str} > {sam_file} 2> {bwa_log}"
                )
                logging.info(f"[{sample_id}] BWA MEM (threads={threads})")
                if not run_cmd(bwa_cmd) or not sam_file.exists():
                    return
                bam_file = sample_dir / f"{sample_id}.bam"
                sorted_bam = sample_dir / f"{sample_id}.sorted.bam"
                if not run_cmd(f"samtools view -S -b {sam_file} -o {bam_file}"): return
                if not run_cmd(f"samtools sort -@ {threads} {bam_file} -o {sorted_bam}"): return
                run_cmd(f"samtools index {sorted_bam}")
                sorted_bam = str(sorted_bam)


        # ----- Variant calling: create TWO vcfs -----
        if sorted_bam and (not vcf_lineage or not vcf_resist):
            # 1) lineage VCF (targeted): BED on NC_000962.3
            vcf_lineage = str(sample_dir / f"{sample_id}.lineage.vcf")
            tmp_bed = sample_dir / "targets.nc.bed"
            rewrite_bed_contig(bed_file, tmp_bed, target="NC_000962.3")
            if not run_cmd(
                f"gatk HaplotypeCaller "
                f"-R {ref_genome} -I {sorted_bam} -O {vcf_lineage} "
                f"-L {tmp_bed} --sample-ploidy 1"
            ):
                return


            # 2) resistance VCF (whole genome)
            vcf_resist = str(sample_dir / f"{sample_id}.resist.vcf")
            if not run_cmd(
                f"gatk HaplotypeCaller "
                f"-R {ref_genome} -I {sorted_bam} -O {vcf_resist} "
                f"--sample-ploidy 1"
            ):
                return


    # ---------------- lineage scoring (use lineage VCF) ----------------
    sample_snps = set()
    with open(vcf_lineage, 'r') as vcf:
        for line in vcf:
            if line.startswith('#'): continue
            chrom, pos = line.strip().split('\t')[:2]
            chrom = contig_alias_to_nc(chrom)
            try:
                sample_snps.add((chrom, int(pos)))
            except ValueError:
                pass


    results, matched_counts = {}, {}
    for lineage, ref_snps in lineage_references.items():
        matched = len(ref_snps & sample_snps)
        total = len(ref_snps)
        results[lineage] = pct(matched, total)
        matched_counts[lineage] = f"{matched}/{total}"


    formatted = {lin: f"{results.get(lin, 0.0):.2f}%" for lin in lineage_references}
    sorted_lineages = sorted(results.items(), key=lambda x: x[1], reverse=True)


    if not sorted_lineages or sorted_lineages[0][1] == 0.0:
        predicted = "Unpredictable: insufficient SNP evidence"
    else:
        top = sorted_lineages[0][1]
        predicted_set = [lin for lin, pc in sorted_lineages if pc >= 90.0 or pc >= top - 25.0]
        predicted = ";".join(sorted(set(predicted_set)))


    high_conf = [lin for lin, pc in results.items() if pc >= 90.0]
    mixed = ";".join(sorted(high_conf)) if len(high_conf) > 1 else "-"


    elapsed_minutes = (time.time() - start_time) / 60
    lin_out = sample_dir / f"{sample_id}_lineage_result.txt"
    with open(lin_out, 'w') as f:
        f.write(f"Sample              : {sample_id}\n")
        f.write(f"Predicted lineage   : {predicted}\n")
        f.write(f"Mixed lineage       : {mixed}\n\n")
        f.write(f"{'Lineage':<20}{'Probability':<15}{'Matched SNPs'}\n")
        f.write("-" * 50 + "\n")
        for lineage in sorted(lineage_references.keys()):
            f.write(f"{lineage:<20}{formatted.get(lineage, '0.00%'):<15}{matched_counts.get(lineage, '0/0')}\n")
        f.write("\n" + "=" * 60 + "\n")
        f.write(f"Total Runtime: {elapsed_minutes:.2f} minutes\n")
    logging.info(f"[{sample_id}] Lineage result saved to {lin_out}")


    # ---------------- mutation summary (snpEff + WHO) ----------
    if mutations_csv:
        ann_vcf = sample_dir / f"{sample_id}_ann.vcf"
        ann_csv = sample_dir / f"{sample_id}_annotated.csv"
        all_csv = sample_dir / "matched_all.csv"
        res_csv = sample_dir / "matched_mutations.csv"
        res_summary_csv = sample_dir / "mutations_summary.csv"
        final_txt = sample_dir / f"{sample_id}_final_prediction.txt"
        final_csv = sample_dir / f"{sample_id}_summary.csv"


        # Rename contig to 'Chromosome' for snpEff (use RESISTANCE VCF)
        vcf_chr = sample_dir / f"{sample_id}.resist.chr.vcf"
        chrmap = sample_dir / "chr_rename.map"
        write_chrmap(str(chrmap), to_chr=True)  # NC_000962.3 -> Chromosome
        if not run_cmd(f"bcftools annotate --rename-chrs {chrmap} {vcf_resist} -Ov -o {vcf_chr}"):
            logging.error(f"[{sample_id}] Failed to rename contig to 'Chromosome' for snpEff.")
            return


        # snpEff annotate
        if not run_cmd(f'{snpeff_cmd} {snpeff_db} {vcf_chr} > {ann_vcf}'):
            logging.error(f"[{sample_id}] snpEff failed; skipping mutation summary.")
            return


        df = parse_snpeff_vcf_to_df(str(ann_vcf))
        df.to_csv(ann_csv, index=False)
        logging.info(f"[{sample_id}] wrote {ann_csv}")


        match_catalogue(str(ann_csv), mutations_csv, str(all_csv), str(res_csv))


        try:
            _df_dbg = pd.read_csv(res_csv)
            logging.info(f"[{sample_id}] matched_mutations rows: {_df_dbg.shape[0]}  path: {res_csv}")
            if "drug" in _df_dbg.columns:
                logging.info(f"[{sample_id}] drugs seen: {sorted(set(str(x).lower() for x in _df_dbg['drug']))}")
        except Exception as e:
            logging.warning(f"[{sample_id}] could not preview matched_mutations: {e}")


        write_drug_summary(str(res_csv), str(res_summary_csv))
        predict_final_category(str(res_csv), str(final_txt), str(final_csv))


# ------------------------------- CLI ---------------------------
def main():
    ap = argparse.ArgumentParser(
        description="LineageXpress: lineage + mutation summary + WHO category (contig alias aware)",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    ap.add_argument("--sample_list", "--fastq", dest="sample_list", required=True,
                    help="Text file: one sample spec per line ('R1,R2' | file/dir path | bare ID)")
    ap.add_argument("--ref_genome", default=DEFAULT_REF, help="Reference FASTA")
    ap.add_argument("--output_dir", default="results", help="Output directory")
    ap.add_argument("--snp_file", default=DEFAULT_SNP_FILE, help="Lineage SNP reference file")
    ap.add_argument("--bed_file", default=DEFAULT_BED_FILE, help="Target regions BED file for lineage VCF")
    ap.add_argument("--bam_file", default=None, help="Override: use this BAM for all samples")
    ap.add_argument("--vcf_file", default=None, help="Override: use this VCF for all samples")
    ap.add_argument("--threads", type=int, default=1, help="Threads for BWA / samtools")


    # Mutation summary options
    ap.add_argument("--mutations_csv", required=False,
                    help="TB-Profiler WHO mutations.csv; enables mutation summary if provided")
    ap.add_argument("--snpeff_cmd", default="snpEff",
                    help='Command for snpEff (e.g., "snpEff" or "java -Xmx4g -jar /path/snpEff.jar")')
    ap.add_argument("--snpeff_db", default=DEFAULT_SNPEFF_DB,
                    help="snpEff DB name")


    args = ap.parse_args()


    # Preflight
    for path, label in [(args.ref_genome, "ref_genome"),
                        (args.snp_file, "snp_file"),
                        (args.bed_file, "bed_file")]:
        if path and not os.path.exists(path):
            raise FileNotFoundError(f"{label} not found at '{path}'")


    lineage_refs = load_lineage_snp_file(args.snp_file)
    lineage_refs = lineage_snps_to_target(lineage_refs, target="NC_000962.3")


    if not os.path.exists(args.sample_list):
        raise FileNotFoundError(f"sample_list file not found at '{args.sample_list}'")
    with open(args.sample_list, "r") as f:
        samples = [line.strip() for line in f if line.strip()]


    os.makedirs(args.output_dir, exist_ok=True)


    for sample in samples:
        process_sample(
            sample_path=sample,
            ref_genome=args.ref_genome,
            output_dir=args.output_dir,
            lineage_references=lineage_refs,
            bed_file=args.bed_file,
            threads=args.threads,
            bam_override=args.bam_file,
            vcf_override=args.vcf_file,
            snpeff_cmd=args.snpeff_cmd,
            snpeff_db=args.snpeff_db,
            mutations_csv=args.mutations_csv,
        )


if __name__ == "__main__":
    main()
