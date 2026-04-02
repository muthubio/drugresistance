#!/usr/bin/env python3
import os
import re
import time
import argparse
import logging
import subprocess
from pathlib import Path
from collections import defaultdict

logging.basicConfig(
    format="%(asctime)s - %(levelname)s - %(message)s",
    level=logging.INFO
)

BASE_DIR = Path(__file__).resolve().parent
DATA_DIR = BASE_DIR / "data"
DEFAULT_REF = str(DATA_DIR / "h37rv.fa")
DEFAULT_SNP_FILE = str(DATA_DIR / "lineage_snp_updated_au13.tsv")
SAMPLE_INPUT_DIR = BASE_DIR / "sample_input"
SAMPLE_DATA_DIR = BASE_DIR / "sample_data"

CONTIG_ALIASES = {"Chromosome", "NC_000962.3", "H37Rv", "chrH37Rv"}


def contig_alias_to_nc(chrom: str) -> str:
    return "NC_000962.3" if chrom in CONTIG_ALIASES else chrom


def lineage_snps_to_target(lineage_snps: dict, target: str = "NC_000962.3"):
    converted = {}
    for lin, sset in lineage_snps.items():
        converted[lin] = {(target if c in CONTIG_ALIASES else c, pos) for (c, pos) in sset}
    return converted


def run_cmd(cmd: str) -> bool:
    logging.info(f"Running: {cmd}")
    try:
        subprocess.run(cmd, shell=True, check=True)
        return True
    except subprocess.CalledProcessError as e:
        logging.error(f"Command failed: {cmd}\n{e}")
        return False


def _find_in_dir(dirpath: Path, sid: str):
    pats = [
        (f"{sid}_1.fastq.gz",  f"{sid}_2.fastq.gz"),
        (f"{sid}_R1.fastq.gz", f"{sid}_R2.fastq.gz"),
        (f"{sid}.1.fastq.gz",  f"{sid}.2.fastq.gz"),
        (f"{sid}_1.fastq",     f"{sid}_2.fastq"),
        (f"{sid}_R1.fastq",    f"{sid}_R2.fastq"),
        (f"{sid}.1.fastq",     f"{sid}.2.fastq"),
        (f"{sid}_out1.fq.gz",      f"{sid}_out2.fq.gz"),
        (f"{sid}_out1.fastq.gz",   f"{sid}_out2.fastq.gz"),
        (f"{sid}_out1.fq",         f"{sid}_out2.fq"),
        (f"{sid}_out1.fastq",      f"{sid}_out2.fastq"),
        (f"{sid}_host_removed_R1.fastq.gz", f"{sid}_host_removed_R2.fastq.gz"),
        (f"{sid}_host_removed_R1.fq.gz",    f"{sid}_host_removed_R2.fq.gz"),
        (f"{sid}_host_removed_R1.fastq",    f"{sid}_host_removed_R2.fastq"),
        (f"{sid}_host_removed_R1.fq",       f"{sid}_host_removed_R2.fq"),
        (f"{sid}_host_removed_1.fastq.gz",  f"{sid}_host_removed_2.fastq.gz"),
        (f"{sid}_host_removed_1.fq.gz",     f"{sid}_host_removed_2.fq.gz"),
        (f"{sid}_host_removed_1.fastq",     f"{sid}_host_removed_2.fastq"),
        (f"{sid}_host_removed_1.fq",        f"{sid}_host_removed_2.fq"),
    ]
    for r1, r2 in pats:
        r1p, r2p = dirpath / r1, dirpath / r2
        if r1p.exists() and r2p.exists():
            return "paired", [str(r1p), str(r2p)]
    for ext in (".fastq.gz", ".fq.gz", ".fastq", ".fq"):
        for se in [dirpath / f"{sid}{ext}", dirpath / f"{sid}_host_removed{ext}"]:
            if se.exists():
                return "single", [str(se)]
    try:
        fqs = sorted([
            p for p in dirpath.iterdir()
            if p.is_file() and (p.suffix in {".fastq", ".fq"} or p.name.endswith(".fastq.gz") or p.name.endswith(".fq.gz"))
        ])
    except FileNotFoundError:
        fqs = []
    if len(fqs) == 2:
        return "paired", [str(fqs[0]), str(fqs[1])]
    if len(fqs) == 1:
        return "single", [str(fqs[0])]
    return "unknown", []


def resolve_inputs(sample_spec: str):
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
        if low.endswith((".vcf", ".vcf.gz")):
            return "vcf", [str(p)]
        if low.endswith(".bam"):
            return "bam", [str(p)]
        if low.endswith((".fastq", ".fq", ".fastq.gz", ".fq.gz")):
            return "single", [str(p)]
        return "unknown", []
    if p.is_dir():
        kind, files = _find_in_dir(p, p.name)
        if kind != "unknown":
            return kind, files
        return _find_in_dir(p, p.name)
    sid = s
    for base in [Path.cwd(), SAMPLE_INPUT_DIR, SAMPLE_DATA_DIR]:
        kind, files = _find_in_dir(base, sid)
        if kind != "unknown":
            return kind, files
    return "unknown", []


def load_lineage_snp_file(snp_file):
    lineage_snps = defaultdict(set)
    with open(snp_file, "r") as f:
        for line in f:
            if line.strip() and not line.startswith("#"):
                parts = line.strip().split("\t")
                if len(parts) == 3:
                    lineage, rvlocus, pos = parts
                elif len(parts) == 2:
                    lineage, pos = parts
                    rvlocus = "NC_000962.3"
                else:
                    continue
                try:
                    lineage_snps[lineage].add((rvlocus, int(pos)))
                except ValueError:
                    pass
    logging.info(f"Loaded SNP reference data for {len(lineage_snps)} lineages from {snp_file}")
    return lineage_snps


def pct(matched, total):
    return round((matched / total) * 100, 2) if total > 0 else 0.0


def normalize_sample_id(sample_path: str) -> str:
    name = Path(sample_path).name
    for suffix in [".fastq.gz", ".fq.gz", ".fastq", ".fq", ".bam", ".vcf.gz", ".vcf"]:
        if name.endswith(suffix):
            name = name[:-len(suffix)]
            break
    name = re.sub(r"_host_removed$", "", name)
    name = re.sub(r"(_host_removed_R1|_host_removed_R2|_host_removed_1|_host_removed_2)$", "", name)
    name = re.sub(r"(_R?1|_R?2|\.R1|\.R2|_1|_2|\.1|\.2|_out1|_out2)$", "", name)
    return name


def parse_vcf_positions(vcf_path: str):
    sample_snps = set()
    opener = open
    mode = "r"
    if str(vcf_path).endswith(".gz"):
        import gzip
        opener = gzip.open
        mode = "rt"
    with opener(vcf_path, mode) as vcf:
        for line in vcf:
            if line.startswith("#"):
                continue
            parts = line.strip().split("\t")
            if len(parts) < 8:
                continue
            chrom, pos, _id, ref, alt, _qual, _filter, info = parts[:8]
            chrom = contig_alias_to_nc(chrom)
            if len(ref) == 1 and all(len(a) == 1 for a in alt.split(",")):
                dp = None
                mq = None
                for item in info.split(";"):
                    if item.startswith("DP="):
                        try:
                            dp = int(item.split("=", 1)[1])
                        except ValueError:
                            pass
                    elif item.startswith("MQ="):
                        try:
                            mq = round(float(item.split("=", 1)[1]), 2)
                        except ValueError:
                            pass
                try:
                    sample_snps.add((chrom, int(pos), dp, mq))
                except ValueError:
                    pass
    return sample_snps


def evaluate_lineages(lineage_references, sample_snps):
    sample_pos_map = {(chrom, pos): {"dp": dp, "mq": mq} for chrom, pos, dp, mq in sample_snps}
    sample_pos_set = set(sample_pos_map.keys())

    results = {}
    matched_counts = {}
    matched_positions = {}

    for lineage, ref_snps in lineage_references.items():
        lineage_positions = set(ref_snps)
        matched = sorted(lineage_positions & sample_pos_set, key=lambda x: x[1])
        total = len(lineage_positions)
        results[lineage] = pct(len(matched), total)
        matched_counts[lineage] = f"{len(matched)}/{total}"
        matched_positions[lineage] = matched

    return results, matched_counts, matched_positions, sample_pos_map


def choose_prediction(results, matched_positions):
    lineages_with_hits = [(lin, len(pos_list), results.get(lin, 0.0)) for lin, pos_list in matched_positions.items() if pos_list]
    if not lineages_with_hits:
        return "Unpredictable: insufficient SNP evidence", []

    lineages_with_hits.sort(key=lambda x: (x[1], x[2], x[0]), reverse=True)
    best_count = lineages_with_hits[0][1]
    best_prob = lineages_with_hits[0][2]
    best_lineages = [lin for lin, cnt, prob in lineages_with_hits if cnt == best_count and prob == best_prob]

    if best_count == 1:
        if len(best_lineages) == 1:
            return best_lineages[0], best_lineages
        if len(best_lineages) <= 3:
            return ";".join(sorted(best_lineages)), best_lineages
        return "Mixed lineage: ambiguous SNP", best_lineages

    if len(best_lineages) == 1:
        return best_lineages[0], best_lineages
    if len(best_lineages) <= 3:
        return ";".join(sorted(best_lineages)), best_lineages
    return "Mixed lineage: ambiguous SNP", best_lineages


def summarize_best_snp(best_lineages, matched_positions, sample_pos_map):
    if not best_lineages:
        return "NA", "NA", "NA", "NA"

    pos_counter = defaultdict(int)
    for lin in best_lineages:
        for chrom, pos in matched_positions.get(lin, []):
            pos_counter[(chrom, pos)] += 1

    if not pos_counter:
        return "NA", "NA", "NA", "NA"

    best_pos = sorted(pos_counter.items(), key=lambda x: (x[1], x[0][1]), reverse=True)[0][0]
    chrom, pos = best_pos
    metrics = sample_pos_map.get((chrom, pos), {})
    dp = metrics.get("dp", "NA")
    mq = metrics.get("mq", "NA")
    return str(pos), "NA", str(dp) if dp is not None else "NA", str(mq) if mq is not None else "NA"


def write_lineage_report(sample_id, lineage_references, results, matched_counts, matched_positions, sample_pos_map, out_path, elapsed_minutes):
    formatted = {lin: f"{results.get(lin, 0.0):.2f}%" for lin in lineage_references}
    predicted, best_lineages = choose_prediction(results, matched_positions)
    high_conf = [lin for lin, pc in results.items() if pc >= 90.0]
    mixed = ";".join(sorted(high_conf)) if len(high_conf) > 1 else "-"
    best_pos, base_qual, depth, mapping_quality = summarize_best_snp(best_lineages, matched_positions, sample_pos_map)

    with open(out_path, "w") as f:
        f.write(f"Sample              : {sample_id}\n")
        f.write(f"Predicted lineage   : {predicted}\n")
        f.write(f"Mixed lineage       : {mixed}\n")
        f.write(f"Best SNP position   : {best_pos}\n")
        f.write(f"Base quality        : {base_qual}\n")
        f.write(f"Depth               : {depth}\n")
        f.write(f"Mapping quality     : {mapping_quality}\n\n")
        f.write(f"{'Lineage':<20}{'Probability':<15}{'Matched SNPs':<15}{'Matched Positions'}\n")
        f.write("-" * 90 + "\n")
        for lineage in sorted(lineage_references.keys()):
            pos_str = ",".join(str(pos) for _chrom, pos in matched_positions.get(lineage, [])) if matched_positions.get(lineage) else "-"
            f.write(f"{lineage:<20}{formatted.get(lineage, '0.00%'):<15}{matched_counts.get(lineage, '0/0'):<15}{pos_str}\n")
        f.write("\n" + "=" * 60 + "\n")
        f.write(f"Total Runtime: {elapsed_minutes:.2f} minutes\n")


def call_variants_bcftools(sorted_bam: str, ref_genome: str, out_vcf_gz: str, threads: int = 1, min_bq: int = 20,
                           min_mq: int = 30, max_depth: int = 100000):
    cmd = (
        f"bcftools mpileup --threads {threads} -f {ref_genome} -q {min_mq} -Q {min_bq} -d {max_depth} -Ou {sorted_bam} | "
        f"bcftools call --threads {threads} -mv --ploidy 1 -Oz -o {out_vcf_gz}"
    )
    if not run_cmd(cmd):
        return False
    return run_cmd(f"bcftools index -f {out_vcf_gz}")


def process_sample(sample_spec, ref_genome, output_dir, lineage_references, threads=1,
                   bam_override=None, vcf_override=None, min_bq=20, min_mq=30, max_depth=100000):
    start_time = time.time()
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    logs_dir = output_dir / "logs"
    logs_dir.mkdir(exist_ok=True)
    sample_id = normalize_sample_id(sample_spec)
    sample_dir = output_dir / sample_id
    sample_dir.mkdir(exist_ok=True)

    sorted_bam = None
    vcf_lineage = None

    if vcf_override and Path(vcf_override).exists():
        vcf_lineage = str(Path(vcf_override))
        logging.info(f"[{sample_id}] Using provided VCF: {vcf_lineage}")
    elif bam_override and Path(bam_override).exists():
        sorted_bam = str(Path(bam_override))
        logging.info(f"[{sample_id}] Using provided BAM: {sorted_bam}")
    else:
        kind, inputs = resolve_inputs(str(sample_spec))
        if kind == "unknown":
            logging.error(f"[{sample_id}] Could not resolve inputs for '{sample_spec}'")
            return
        if kind == "vcf":
            vcf_lineage = inputs[0]
            logging.info(f"[{sample_id}] Detected VCF: {vcf_lineage}")
        elif kind == "bam":
            sorted_bam = inputs[0]
            logging.info(f"[{sample_id}] Detected BAM: {sorted_bam}")
        else:
            bwa_log = logs_dir / f"{sample_id}.bwa.log"
            sam_file = sample_dir / f"{sample_id}.sam"
            bam_file = sample_dir / f"{sample_id}.bam"
            sorted_bam_path = sample_dir / f"{sample_id}.sorted.bam"
            fastq_str = " ".join(inputs)
            bwa_cmd = (
                f"bwa mem -t {threads} -R '@RG\\tID:{sample_id}\\tSM:{sample_id}\\tPL:ILLUMINA' "
                f"{ref_genome} {fastq_str} > {sam_file} 2> {bwa_log}"
            )
            logging.info(f"[{sample_id}] Mapping host-removed reads to H37Rv")
            if not run_cmd(bwa_cmd):
                return
            if not run_cmd(f"samtools view -@ {threads} -bS {sam_file} -o {bam_file}"):
                return
            if not run_cmd(f"samtools sort -@ {threads} {bam_file} -o {sorted_bam_path}"):
                return
            if not run_cmd(f"samtools index {sorted_bam_path}"):
                return
            if Path(sam_file).exists():
                Path(sam_file).unlink()
            if Path(bam_file).exists():
                Path(bam_file).unlink()
            sorted_bam = str(sorted_bam_path)

    if sorted_bam and not vcf_lineage:
        vcf_lineage = str(sample_dir / f"{sample_id}.bcftools.vcf.gz")
        logging.info(f"[{sample_id}] Calling variants with bcftools mpileup/call")
        if not call_variants_bcftools(sorted_bam, ref_genome, vcf_lineage, threads, min_bq, min_mq, max_depth):
            return

    sample_snps = parse_vcf_positions(vcf_lineage)
    results, matched_counts, matched_positions, sample_pos_map = evaluate_lineages(lineage_references, sample_snps)

    elapsed_minutes = (time.time() - start_time) / 60
    result_txt = sample_dir / f"{sample_id}_lineage_result.txt"
    write_lineage_report(sample_id, lineage_references, results, matched_counts, matched_positions, sample_pos_map, result_txt, elapsed_minutes)
    logging.info(f"[{sample_id}] Lineage result saved to {result_txt}")


def main():
    ap = argparse.ArgumentParser(
        description="Metagenomics lineage pipeline with matched SNP reporting and controlled mixed-lineage output",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    ap.add_argument("--sample_list", "--fastq", dest="sample_list", required=True,
                    help="Text file: one sample spec per line ('R1,R2' | FASTQ/BAM/VCF path | directory | bare ID)")
    ap.add_argument("--ref_genome", default=DEFAULT_REF, help="Reference FASTA (H37Rv)")
    ap.add_argument("--output_dir", default="results", help="Output directory")
    ap.add_argument("--snp_file", default=DEFAULT_SNP_FILE, help="LineageXpress SNP reference file")
    ap.add_argument("--bam_file", default=None, help="Override: use this BAM for all samples")
    ap.add_argument("--vcf_file", default=None, help="Override: use this VCF for all samples")
    ap.add_argument("--threads", type=int, default=1, help="Threads for bwa/samtools/bcftools")
    ap.add_argument("--min_bq", type=int, default=20, help="Minimum base quality for mpileup")
    ap.add_argument("--min_mq", type=int, default=30, help="Minimum mapping quality for mpileup")
    ap.add_argument("--max_depth", type=int, default=100000, help="Maximum read depth for mpileup")
    args = ap.parse_args()

    for path, label in [(args.ref_genome, "ref_genome"), (args.snp_file, "snp_file")]:
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
        process_sample(sample, args.ref_genome, args.output_dir, lineage_refs,
                       args.threads, args.bam_file, args.vcf_file,
                       args.min_bq, args.min_mq, args.max_depth)


if __name__ == "__main__":
    main()
