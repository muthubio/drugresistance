#!/usr/bin/env python3
"""
LineageXpress v13
==================

Fast MTBC lineage and drug-resistance prediction pipeline with:
  1. sample-level parallel execution
  2. targeted BED variant calling for lineage prediction
  3. targeted BED variant calling for resistance prediction
  4. FASTQ, BAM, and VCF input support
  5. bcftools-based targeted variant calling by default

Recommended strategy:
  - Map reads to the full H37Rv reference.
  - Restrict variant calling to targeted BED regions.
  - Use lineage BED and resistance BED separately, or provide one combined BED.

Why this design?
  Mapping to the full H37Rv genome keeps original coordinates stable.
  BED-restricted variant calling avoids scanning/calling the complete genome.
  This is faster and safer than mapping directly to small gene FASTA files.

Example sample sheet TSV:
  sample_id\tinput_type\tr1\tr2\tbam\tvcf
  S1\tfastq\t/path/S1_R1.fastq.gz\t/path/S1_R2.fastq.gz\t\t
  S2\tbam\t\t\t/path/S2.sorted.bam\t
  S3\tvcf\t\t\t\t/path/S3.vcf.gz

Example run:
  python lineageexpress_new_version13.py \
    --sample-sheet samples.tsv \
    --ref /home/muthukumarb/drugresistance/data/reference/h37rv.fa \
    --lineage-db /home/muthukumarb/drugresistance/data/database/lineage_snp_updated_au13.tsv \
    --resistance-db /home/muthukumarb/drugresistance/data/database/resistance_mutations.tsv \
    --lineage-bed /home/muthukumarb/drugresistance/data/database/lineage_targets.bed \
    --resistance-bed /home/muthukumarb/drugresistance/data/database/resistance_targets.bed \
    --outdir results_v13 \
    --max-samples-parallel 4 \
    --threads-per-sample 8

Database expectations:
  Lineage DB should contain at least:
    position, ref, alt, lineage
  Optional columns are allowed.

  Resistance DB should contain at least:
    position, ref, alt, drug
  Optional columns such as gene, mutation, confidence are allowed.

The script tries to auto-detect common column names.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import os
import shutil
import subprocess
import sys
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Tuple


# -----------------------------
# Data structures
# -----------------------------

@dataclass
class Sample:
    sample_id: str
    input_type: str
    r1: str = ""
    r2: str = ""
    bam: str = ""
    vcf: str = ""


@dataclass
class Config:
    ref: str
    lineage_db: str
    resistance_db: str
    lineage_bed: str
    resistance_bed: str
    outdir: str
    threads_per_sample: int
    caller: str
    min_dp: int
    min_qual: float
    keep_intermediate: bool
    bwa: str
    samtools: str
    bcftools: str
    gatk: str


# -----------------------------
# Utilities
# -----------------------------

def log(message: str) -> None:
    ts = time.strftime("%Y-%m-%d %H:%M:%S")
    print(f"[{ts}] [INFO] {message}", flush=True)


def fail(message: str) -> None:
    ts = time.strftime("%Y-%m-%d %H:%M:%S")
    print(f"[{ts}] [ERROR] {message}", file=sys.stderr, flush=True)
    raise RuntimeError(message)


def run_cmd(cmd: List[str], log_file: Path, cwd: Optional[Path] = None) -> None:
    log_file.parent.mkdir(parents=True, exist_ok=True)
    with log_file.open("a") as lf:
        lf.write("\n[CMD] " + " ".join(cmd) + "\n")
        lf.flush()
        p = subprocess.run(cmd, stdout=lf, stderr=lf, cwd=str(cwd) if cwd else None)
    if p.returncode != 0:
        raise RuntimeError(f"Command failed with exit code {p.returncode}: {' '.join(cmd)}. See log: {log_file}")


def check_executable(exe: str) -> None:
    if shutil.which(exe) is None:
        fail(f"Required executable not found in PATH: {exe}")


def check_file(path: str, label: str) -> None:
    if path and not Path(path).exists():
        fail(f"{label} file not found: {path}")


def open_text_auto(path: str):
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "r")


def detect_delimiter(path: str) -> str:
    with open_text_auto(path) as f:
        first = f.readline()
    if "\t" in first:
        return "\t"
    return ","


def normalize_header(name: str) -> str:
    return name.strip().lower().replace(" ", "_").replace("-", "_")


def get_first(row: Dict[str, str], candidates: Iterable[str], default: str = "") -> str:
    for c in candidates:
        if c in row and row[c] not in (None, ""):
            return str(row[c]).strip()
    return default


# -----------------------------
# Input parsing
# -----------------------------

def read_sample_sheet(path: str) -> List[Sample]:
    check_file(path, "sample sheet")
    delim = detect_delimiter(path)
    samples: List[Sample] = []
    with open_text_auto(path) as f:
        reader = csv.DictReader(f, delimiter=delim)
        if reader.fieldnames is None:
            fail("Sample sheet has no header")
        reader.fieldnames = [normalize_header(x) for x in reader.fieldnames]
        for row in reader:
            row = {normalize_header(k): (v.strip() if isinstance(v, str) else v) for k, v in row.items()}
            sample_id = get_first(row, ["sample_id", "sample", "id", "name"])
            input_type = get_first(row, ["input_type", "type", "mode"]).lower()
            if not sample_id:
                fail("Sample sheet row missing sample_id")
            if input_type not in {"fastq", "bam", "vcf"}:
                fail(f"Sample {sample_id}: input_type must be fastq, bam, or vcf; got '{input_type}'")
            samples.append(
                Sample(
                    sample_id=sample_id,
                    input_type=input_type,
                    r1=get_first(row, ["r1", "read1", "fastq1", "fq1"]),
                    r2=get_first(row, ["r2", "read2", "fastq2", "fq2"]),
                    bam=get_first(row, ["bam", "bam_path"]),
                    vcf=get_first(row, ["vcf", "vcf_path"]),
                )
            )
    if not samples:
        fail("No samples found in sample sheet")
    return samples


def validate_sample(sample: Sample) -> None:
    if sample.input_type == "fastq":
        check_file(sample.r1, f"{sample.sample_id} R1 FASTQ")
        if sample.r2:
            check_file(sample.r2, f"{sample.sample_id} R2 FASTQ")
    elif sample.input_type == "bam":
        check_file(sample.bam, f"{sample.sample_id} BAM")
    elif sample.input_type == "vcf":
        check_file(sample.vcf, f"{sample.sample_id} VCF")


# -----------------------------
# BED handling
# -----------------------------

def combine_beds(lineage_bed: str, resistance_bed: str, out_bed: Path) -> str:
    """Create a combined BED file from lineage and resistance BEDs."""
    out_bed.parent.mkdir(parents=True, exist_ok=True)
    seen = set()
    with out_bed.open("w") as out:
        for bed in [lineage_bed, resistance_bed]:
            if not bed:
                continue
            check_file(bed, "BED")
            with open_text_auto(bed) as f:
                for line in f:
                    line = line.strip()
                    if not line or line.startswith("#"):
                        continue
                    if line not in seen:
                        seen.add(line)
                        out.write(line + "\n")
    return str(out_bed)


# -----------------------------
# Mapping and variant calling
# -----------------------------

def ensure_bam_index(bam: str, cfg: Config, log_file: Path) -> None:
    bai1 = Path(bam + ".bai")
    bai2 = Path(str(Path(bam).with_suffix("")) + ".bai")
    if not bai1.exists() and not bai2.exists():
        run_cmd([cfg.samtools, "index", "-@", str(cfg.threads_per_sample), bam], log_file)


def map_fastq_to_bam(sample: Sample, sample_dir: Path, cfg: Config, log_file: Path) -> str:
    bam_dir = sample_dir / "bam"
    bam_dir.mkdir(parents=True, exist_ok=True)
    sam = bam_dir / f"{sample.sample_id}.sam"
    sorted_bam = bam_dir / f"{sample.sample_id}.sorted.bam"

    if sorted_bam.exists() and Path(str(sorted_bam) + ".bai").exists():
        log(f"{sample.sample_id}: sorted BAM already exists; skipping mapping")
        return str(sorted_bam)

    if sample.r2:
        bwa_cmd = [cfg.bwa, "mem", "-t", str(cfg.threads_per_sample), cfg.ref, sample.r1, sample.r2]
    else:
        bwa_cmd = [cfg.bwa, "mem", "-t", str(cfg.threads_per_sample), cfg.ref, sample.r1]

    sort_cmd = [cfg.samtools, "sort", "-@", str(cfg.threads_per_sample), "-o", str(sorted_bam), "-"]

    with log_file.open("a") as lf:
        lf.write("\n[CMD] " + " ".join(bwa_cmd) + " | " + " ".join(sort_cmd) + "\n")
        lf.flush()
        p1 = subprocess.Popen(bwa_cmd, stdout=subprocess.PIPE, stderr=lf)
        p2 = subprocess.Popen(sort_cmd, stdin=p1.stdout, stdout=lf, stderr=lf)
        if p1.stdout:
            p1.stdout.close()
        rc2 = p2.wait()
        rc1 = p1.wait()
    if rc1 != 0 or rc2 != 0:
        raise RuntimeError(f"{sample.sample_id}: BWA/Samtools mapping failed. See log: {log_file}")

    ensure_bam_index(str(sorted_bam), cfg, log_file)
    return str(sorted_bam)


def call_variants_bcftools(sample_id: str, bam: str, bed: str, out_vcf: Path, cfg: Config, log_file: Path) -> str:
    out_vcf.parent.mkdir(parents=True, exist_ok=True)
    if out_vcf.exists() and Path(str(out_vcf) + ".tbi").exists():
        log(f"{sample_id}: VCF already exists; skipping variant calling: {out_vcf}")
        return str(out_vcf)

    mpileup_cmd = [
        cfg.bcftools, "mpileup",
        "--threads", str(cfg.threads_per_sample),
        "-f", cfg.ref,
        "-R", bed,
        "-Ou",
        bam,
    ]
    call_cmd = [
        cfg.bcftools, "call",
        "--threads", str(cfg.threads_per_sample),
        "-mv",
        "-Oz",
        "-o", str(out_vcf),
    ]

    with log_file.open("a") as lf:
        lf.write("\n[CMD] " + " ".join(mpileup_cmd) + " | " + " ".join(call_cmd) + "\n")
        lf.flush()
        p1 = subprocess.Popen(mpileup_cmd, stdout=subprocess.PIPE, stderr=lf)
        p2 = subprocess.Popen(call_cmd, stdin=p1.stdout, stdout=lf, stderr=lf)
        if p1.stdout:
            p1.stdout.close()
        rc2 = p2.wait()
        rc1 = p1.wait()
    if rc1 != 0 or rc2 != 0:
        raise RuntimeError(f"{sample_id}: bcftools variant calling failed. See log: {log_file}")

    run_cmd([cfg.bcftools, "index", "-t", "-f", str(out_vcf)], log_file)
    return str(out_vcf)


def call_variants_gatk(sample_id: str, bam: str, bed: str, out_vcf: Path, cfg: Config, log_file: Path) -> str:
    out_vcf.parent.mkdir(parents=True, exist_ok=True)
    if out_vcf.exists() and Path(str(out_vcf) + ".tbi").exists():
        log(f"{sample_id}: VCF already exists; skipping GATK calling: {out_vcf}")
        return str(out_vcf)

    cmd = [
        cfg.gatk, "HaplotypeCaller",
        "-R", cfg.ref,
        "-I", bam,
        "-L", bed,
        "-O", str(out_vcf),
    ]
    run_cmd(cmd, log_file)
    return str(out_vcf)


# -----------------------------
# VCF parsing
# -----------------------------

def parse_info_dp(info: str) -> Optional[int]:
    for item in info.split(";"):
        if item.startswith("DP="):
            try:
                return int(float(item.split("=", 1)[1]))
            except ValueError:
                return None
    return None


def load_vcf_variants(vcf_path: str, min_dp: int, min_qual: float) -> Dict[Tuple[int, str, str], Dict[str, str]]:
    variants: Dict[Tuple[int, str, str], Dict[str, str]] = {}
    with open_text_auto(vcf_path) as f:
        for line in f:
            if not line or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 8:
                continue
            chrom, pos, vid, ref, alt, qual, flt, info = parts[:8]
            try:
                pos_i = int(pos)
            except ValueError:
                continue
            try:
                qual_f = float(qual) if qual not in {".", ""} else 0.0
            except ValueError:
                qual_f = 0.0
            if qual_f < min_qual:
                continue
            dp = parse_info_dp(info)
            if dp is not None and dp < min_dp:
                continue
            for alt_allele in alt.split(","):
                key = (pos_i, ref.upper(), alt_allele.upper())
                variants[key] = {
                    "chrom": chrom,
                    "pos": str(pos_i),
                    "id": vid,
                    "ref": ref,
                    "alt": alt_allele,
                    "qual": str(qual_f),
                    "filter": flt,
                    "info": info,
                    "dp": str(dp) if dp is not None else "",
                }
    return variants


# -----------------------------
# DB parsing and prediction
# -----------------------------

def load_marker_db(path: str, db_type: str) -> List[Dict[str, str]]:
    if not path:
        return []
    check_file(path, f"{db_type} database")
    delim = detect_delimiter(path)
    rows: List[Dict[str, str]] = []
    with open_text_auto(path) as f:
        reader = csv.DictReader(f, delimiter=delim)
        if reader.fieldnames is None:
            fail(f"{db_type} database has no header")
        reader.fieldnames = [normalize_header(x) for x in reader.fieldnames]
        for row in reader:
            norm = {normalize_header(k): (v.strip() if isinstance(v, str) else "") for k, v in row.items()}
            rows.append(norm)
    return rows


def marker_key(row: Dict[str, str]) -> Optional[Tuple[int, str, str]]:
    pos_s = get_first(row, ["position", "pos", "genome_position", "snp_position", "coordinate"])
    ref = get_first(row, ["ref", "reference", "reference_allele", "ref_allele"])
    alt = get_first(row, ["alt", "alternate", "alternate_allele", "alt_allele", "mut", "mutation_allele"])
    if not pos_s or not ref or not alt:
        return None
    try:
        pos = int(float(pos_s))
    except ValueError:
        return None
    return (pos, ref.upper(), alt.upper())


def predict_lineage(variants: Dict[Tuple[int, str, str], Dict[str, str]], lineage_db: List[Dict[str, str]]) -> Tuple[str, List[Dict[str, str]]]:
    hits: List[Dict[str, str]] = []
    lineage_counts: Dict[str, int] = {}
    for row in lineage_db:
        key = marker_key(row)
        if key and key in variants:
            lineage = get_first(row, ["lineage", "sublineage", "lineage_name", "class", "label"], "Unknown")
            hit = dict(row)
            hit.update({
                "matched_pos": str(key[0]),
                "matched_ref": key[1],
                "matched_alt": key[2],
                "variant_qual": variants[key].get("qual", ""),
                "variant_dp": variants[key].get("dp", ""),
            })
            hits.append(hit)
            lineage_counts[lineage] = lineage_counts.get(lineage, 0) + 1

    if not lineage_counts:
        return "Unclassified", hits

    best = sorted(lineage_counts.items(), key=lambda x: (-x[1], x[0]))[0][0]
    return best, hits


def predict_resistance(variants: Dict[Tuple[int, str, str], Dict[str, str]], resistance_db: List[Dict[str, str]]) -> Tuple[str, List[Dict[str, str]]]:
    hits: List[Dict[str, str]] = []
    drugs = set()
    for row in resistance_db:
        key = marker_key(row)
        if key and key in variants:
            drug = get_first(row, ["drug", "antibiotic", "drug_name"], "Unknown")
            gene = get_first(row, ["gene", "locus", "target_gene"], "")
            mutation = get_first(row, ["mutation", "aa_change", "change", "variant"], "")
            confidence = get_first(row, ["confidence", "grade", "evidence"], "")
            hit = dict(row)
            hit.update({
                "matched_pos": str(key[0]),
                "matched_ref": key[1],
                "matched_alt": key[2],
                "drug": drug,
                "gene": gene,
                "mutation": mutation,
                "confidence": confidence,
                "variant_qual": variants[key].get("qual", ""),
                "variant_dp": variants[key].get("dp", ""),
            })
            hits.append(hit)
            drugs.add(drug)

    if not hits:
        return "No known resistance mutation detected", hits
    return ";".join(sorted(drugs)), hits


# -----------------------------
# Output writers
# -----------------------------

def write_tsv(path: Path, rows: List[Dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        with path.open("w") as f:
            f.write("No_hits\n")
        return
    headers = sorted({k for row in rows for k in row.keys()})
    with path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=headers, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def write_summary(path: Path, summary_rows: List[Dict[str, str]]) -> None:
    headers = [
        "sample_id", "input_type", "status", "lineage_prediction", "resistance_prediction",
        "lineage_hits", "resistance_hits", "vcf_used", "runtime_seconds", "error"
    ]
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=headers, delimiter="\t")
        writer.writeheader()
        writer.writerows(summary_rows)


# -----------------------------
# Per-sample pipeline
# -----------------------------

def process_one_sample(sample: Sample, cfg: Config) -> Dict[str, str]:
    start = time.time()
    sample_dir = Path(cfg.outdir) / sample.sample_id
    sample_dir.mkdir(parents=True, exist_ok=True)
    log_file = sample_dir / "logs" / f"{sample.sample_id}.log"

    result = {
        "sample_id": sample.sample_id,
        "input_type": sample.input_type,
        "status": "FAILED",
        "lineage_prediction": "",
        "resistance_prediction": "",
        "lineage_hits": "0",
        "resistance_hits": "0",
        "vcf_used": "",
        "runtime_seconds": "0",
        "error": "",
    }

    try:
        validate_sample(sample)
        log(f"{sample.sample_id}: started")

        combined_bed = combine_beds(
            cfg.lineage_bed,
            cfg.resistance_bed,
            sample_dir / "targets" / f"{sample.sample_id}.combined_targets.bed",
        )

        if sample.input_type == "vcf":
            vcf_path = sample.vcf
        else:
            if sample.input_type == "fastq":
                bam = map_fastq_to_bam(sample, sample_dir, cfg, log_file)
            else:
                bam = sample.bam
                ensure_bam_index(bam, cfg, log_file)

            vcf_path = str(sample_dir / "vcf" / f"{sample.sample_id}.targeted.vcf.gz")
            if cfg.caller == "bcftools":
                call_variants_bcftools(sample.sample_id, bam, combined_bed, Path(vcf_path), cfg, log_file)
            elif cfg.caller == "gatk":
                call_variants_gatk(sample.sample_id, bam, combined_bed, Path(vcf_path), cfg, log_file)
            else:
                fail(f"Unsupported caller: {cfg.caller}")

        variants = load_vcf_variants(vcf_path, cfg.min_dp, cfg.min_qual)
        lineage_db = load_marker_db(cfg.lineage_db, "lineage")
        resistance_db = load_marker_db(cfg.resistance_db, "resistance") if cfg.resistance_db else []

        lineage_prediction, lineage_hits = predict_lineage(variants, lineage_db)
        resistance_prediction, resistance_hits = predict_resistance(variants, resistance_db)

        write_tsv(sample_dir / f"{sample.sample_id}_lineage_hits.tsv", lineage_hits)
        write_tsv(sample_dir / f"{sample.sample_id}_resistance_hits.tsv", resistance_hits)

        with (sample_dir / f"{sample.sample_id}_final_report.txt").open("w") as f:
            f.write(f"Sample_ID\t{sample.sample_id}\n")
            f.write(f"Input_type\t{sample.input_type}\n")
            f.write(f"VCF_used\t{vcf_path}\n")
            f.write(f"Lineage_prediction\t{lineage_prediction}\n")
            f.write(f"Lineage_hits\t{len(lineage_hits)}\n")
            f.write(f"Resistance_prediction\t{resistance_prediction}\n")
            f.write(f"Resistance_hits\t{len(resistance_hits)}\n")
            f.write(f"Mode\tFull-reference mapping + targeted BED variant calling\n")

        result.update({
            "status": "SUCCESS",
            "lineage_prediction": lineage_prediction,
            "resistance_prediction": resistance_prediction,
            "lineage_hits": str(len(lineage_hits)),
            "resistance_hits": str(len(resistance_hits)),
            "vcf_used": vcf_path,
        })
        log(f"{sample.sample_id}: completed")

    except Exception as e:
        result["error"] = str(e)
        with log_file.open("a") as lf:
            lf.write(f"\n[ERROR] {e}\n")
        log(f"{sample.sample_id}: failed: {e}")

    finally:
        result["runtime_seconds"] = str(round(time.time() - start, 2))

    return result


# -----------------------------
# CLI
# -----------------------------

def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="LineageXpress v13: parallel targeted-BED MTBC lineage and resistance prediction"
    )
    p.add_argument("--sample-sheet", required=True, help="TSV/CSV sample sheet with sample_id,input_type,r1,r2,bam,vcf columns")
    p.add_argument("--ref", required=True, help="Full H37Rv reference FASTA")
    p.add_argument("--lineage-db", required=True, help="Lineage SNP database TSV/CSV")
    p.add_argument("--resistance-db", default="", help="Resistance mutation database TSV/CSV")
    p.add_argument("--lineage-bed", required=True, help="BED file covering lineage marker regions")
    p.add_argument("--resistance-bed", default="", help="BED file covering resistance genes/markers")
    p.add_argument("--outdir", required=True, help="Output directory")
    p.add_argument("--max-samples-parallel", type=int, default=1, help="Number of samples to run in parallel")
    p.add_argument("--threads-per-sample", type=int, default=4, help="Threads used inside each sample")
    p.add_argument("--caller", choices=["bcftools", "gatk"], default="bcftools", help="Variant caller for targeted regions")
    p.add_argument("--min-dp", type=int, default=5, help="Minimum INFO/DP for accepting VCF variant")
    p.add_argument("--min-qual", type=float, default=0.0, help="Minimum QUAL for accepting VCF variant")
    p.add_argument("--keep-intermediate", action="store_true", help="Keep intermediate files")
    p.add_argument("--bwa", default="bwa", help="bwa executable")
    p.add_argument("--samtools", default="samtools", help="samtools executable")
    p.add_argument("--bcftools", default="bcftools", help="bcftools executable")
    p.add_argument("--gatk", default="gatk", help="gatk executable")
    return p.parse_args()


def main() -> None:
    args = parse_args()

    os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
    os.environ.setdefault("OMP_NUM_THREADS", "1")
    os.environ.setdefault("MKL_NUM_THREADS", "1")
    os.environ.setdefault("NUMEXPR_NUM_THREADS", "1")

    check_file(args.ref, "reference FASTA")
    check_file(args.lineage_db, "lineage database")
    check_file(args.lineage_bed, "lineage BED")
    if args.resistance_db:
        check_file(args.resistance_db, "resistance database")
    if args.resistance_bed:
        check_file(args.resistance_bed, "resistance BED")

    check_executable(args.samtools)
    if args.caller == "bcftools":
        check_executable(args.bcftools)
    if args.caller == "gatk":
        check_executable(args.gatk)
    check_executable(args.bwa)

    samples = read_sample_sheet(args.sample_sheet)
    Path(args.outdir).mkdir(parents=True, exist_ok=True)

    cfg = Config(
        ref=args.ref,
        lineage_db=args.lineage_db,
        resistance_db=args.resistance_db,
        lineage_bed=args.lineage_bed,
        resistance_bed=args.resistance_bed,
        outdir=args.outdir,
        threads_per_sample=args.threads_per_sample,
        caller=args.caller,
        min_dp=args.min_dp,
        min_qual=args.min_qual,
        keep_intermediate=args.keep_intermediate,
        bwa=args.bwa,
        samtools=args.samtools,
        bcftools=args.bcftools,
        gatk=args.gatk,
    )

    log(f"Loaded {len(samples)} samples")
    log(f"Running up to {args.max_samples_parallel} samples in parallel")
    log(f"Threads per sample: {args.threads_per_sample}")
    log("Mode: full-reference mapping + targeted BED variant calling")

    summary_rows: List[Dict[str, str]] = []
    if args.max_samples_parallel <= 1:
        for sample in samples:
            summary_rows.append(process_one_sample(sample, cfg))
    else:
        with ProcessPoolExecutor(max_workers=args.max_samples_parallel) as ex:
            future_map = {ex.submit(process_one_sample, sample, cfg): sample.sample_id for sample in samples}
            for fut in as_completed(future_map):
                sample_id = future_map[fut]
                try:
                    summary_rows.append(fut.result())
                except Exception as e:
                    summary_rows.append({
                        "sample_id": sample_id,
                        "input_type": "",
                        "status": "FAILED",
                        "lineage_prediction": "",
                        "resistance_prediction": "",
                        "lineage_hits": "0",
                        "resistance_hits": "0",
                        "vcf_used": "",
                        "runtime_seconds": "0",
                        "error": str(e),
                    })

    summary_rows = sorted(summary_rows, key=lambda x: x.get("sample_id", ""))
    write_summary(Path(args.outdir) / "LineageXpress_v13_summary.tsv", summary_rows)

    success = sum(1 for r in summary_rows if r["status"] == "SUCCESS")
    failed = len(summary_rows) - success
    log(f"All done. Success: {success}; Failed: {failed}")
    log(f"Summary: {Path(args.outdir) / 'LineageXpress_v13_summary.tsv'}")

    if failed:
        sys.exit(1)


if __name__ == "__main__":
    main()
