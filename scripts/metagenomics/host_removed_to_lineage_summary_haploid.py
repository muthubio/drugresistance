#!/usr/bin/env python3

import os
import csv
import time
import shlex
import argparse
import logging
import subprocess
from pathlib import Path
from collections import defaultdict
import re
import glob

LOG_FORMAT = "%(asctime)s | %(levelname)s | %(message)s"
logging.basicConfig(level=logging.INFO, format=LOG_FORMAT)
LOGGER = logging.getLogger("host_removed_to_lineage")

FASTQ_EXTENSIONS = [".fastq.gz", ".fq.gz", ".fastq", ".fq"]
PAIR_PATTERNS = [
    ("_R1_001", "_R2_001"), ("_R1", "_R2"), ("_r1", "_r2"), ("_1", "_2"),
    (".R1", ".R2"), (".r1", ".r2"), ("_out1", "_out2"), (".1", ".2"),
    ("_host_removed_R1", "_host_removed_R2"), ("_host_removed_1", "_host_removed_2"),
]
SUMMARY_FIELDS = ["Sample_name", "SNP_position_Lineage_express_db", "Base_quality", "Depth", "Mapping_Quality", "Predicted_lineage"]


def ensure_dir(path):
    Path(path).mkdir(parents=True, exist_ok=True)


def run_command(cmd, shell=False):
    cmd_str = cmd if isinstance(cmd, str) else " ".join(shlex.quote(str(x)) for x in cmd)
    LOGGER.info("CMD: %s", cmd_str)
    proc = subprocess.run(cmd, shell=shell, check=False)
    if proc.returncode != 0:
        raise RuntimeError(f"Command failed ({proc.returncode}): {cmd_str}")


def strip_extension(name):
    for ext in FASTQ_EXTENSIONS:
        if name.endswith(ext):
            return name[:-len(ext)]
    return name


def clean_sample_name(filename):
    base = strip_extension(os.path.basename(filename))
    patterns = [r"_host_removed_R1$", r"_host_removed_R2$", r"_host_removed_1$", r"_host_removed_2$", r"_R1_001$", r"_R2_001$", r"_R1$", r"_R2$", r"_r1$", r"_r2$", r"_1$", r"_2$", r"\.R1$", r"\.R2$", r"\.r1$", r"\.r2$", r"_out1$", r"_out2$", r"\.1$", r"\.2$"]
    for p in patterns:
        base = re.sub(p, "", base)
    return base


def auto_detect_fastq(entry):
    entry = entry.strip()
    if "," in entry:
        parts = [x.strip() for x in entry.split(",") if x.strip()]
        if len(parts) == 2:
            r1, r2 = parts
            return "paired", [r1, r2], clean_sample_name(r1)
    if any(entry.endswith(ext) for ext in FASTQ_EXTENSIONS):
        mate = None
        for left, right in PAIR_PATTERNS:
            if left in entry:
                mate = entry.replace(left, right, 1)
                break
            elif right in entry:
                mate = entry.replace(right, left, 1)
                break
        if mate and os.path.exists(mate):
            if any(x in os.path.basename(entry) for x in ["_R1_001", "_R1", "_r1", "_1", ".R1", ".r1", ".1", "_out1", "_host_removed_R1", "_host_removed_1"]):
                return "paired", [entry, mate], clean_sample_name(entry)
            return "paired", [mate, entry], clean_sample_name(mate)
        return "single", [entry], clean_sample_name(entry)
    parent = os.path.dirname(entry) or "."
    prefix = os.path.basename(entry)
    matches = []
    for ext in FASTQ_EXTENSIONS:
        matches.extend(glob.glob(os.path.join(parent, prefix + "*" + ext)))
    matches = sorted(set(matches))
    if not matches:
        raise FileNotFoundError(f"No FASTQ files found for input: {entry}")
    for f in matches:
        for left, right in PAIR_PATTERNS:
            if left in f:
                mate = f.replace(left, right, 1)
                if mate in matches:
                    return "paired", [f, mate], clean_sample_name(f)
    return "single", [matches[0]], clean_sample_name(matches[0])


def append_file(src, dst):
    if os.path.exists(src) and os.path.getsize(src) > 0:
        with open(src, "r", errors="replace") as s, open(dst, "a") as d:
            for line in s:
                d.write(line)


def get_taxids_from_kraken_report(report_file):
    patterns = ["mycobacterium", "tuberculosis", "africanum", "bovis", "caprae", "canettii", "microti", "pinnipedii", "orygis", "mungi"]
    taxids = set()
    with open(report_file, "r", errors="replace") as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 6:
                continue
            taxid = parts[4].strip()
            name = parts[5].strip().lower()
            if any(pat in name for pat in patterns):
                if taxid:
                    taxids.add(taxid)
    return sorted(taxids)


def fastq_record_iterator(fastq_file):
    with open(fastq_file, "r", errors="replace") as fh:
        while True:
            h = fh.readline()
            if not h:
                break
            s = fh.readline()
            p = fh.readline()
            q = fh.readline()
            if not q:
                break
            yield (h, s, p, q)


def fastq_header_id(header_line):
    h = header_line.strip()
    if h.startswith("@"):
        h = h[1:]
    h = h.split()[0]
    if h.endswith("/1") or h.endswith("/2"):
        h = h[:-2]
    return h


def deduplicate_paired_fastq_by_r1(r1_file, r2_file, out_r1, out_r2):
    seen = set()
    with open(out_r1, "w") as f1_out, open(out_r2, "w") as f2_out:
        for rec1, rec2 in zip(fastq_record_iterator(r1_file), fastq_record_iterator(r2_file)):
            rid = fastq_header_id(rec1[0])
            if rid in seen:
                continue
            seen.add(rid)
            f1_out.writelines(rec1)
            f2_out.writelines(rec2)


def load_lineage_snp_db(snp_db):
    lineage_snps = defaultdict(set)
    with open(snp_db, "r") as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            try:
                if len(parts) == 2:
                    lineage, pos = parts
                    chrom = "NC_000962.3"
                else:
                    lineage, chrom, pos = parts[:3]
                lineage_snps[lineage].add((chrom, int(pos)))
            except Exception:
                continue
    return lineage_snps


def load_bed_positions(bed_file):
    bed_positions = set()
    with open(bed_file, "r") as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            chrom, start, end = line.split("\t")[:3]
            for pos in range(int(start) + 1, int(end) + 1):
                bed_positions.add((chrom, pos))
    return bed_positions


def parse_vcf_variants(vcf_file):
    variant_positions = set()
    stats_by_pos = {}
    opener = ["cat"] if not vcf_file.endswith(".gz") else ["zcat"]
    lines = subprocess.check_output(f"{' '.join(opener)} {shlex.quote(vcf_file)}", shell=True, text=True)
    for line in lines.splitlines():
        if line.startswith("#"):
            continue
        fields = line.split("\t")
        if len(fields) < 8:
            continue
        chrom = fields[0]
        try:
            pos = int(fields[1])
        except ValueError:
            continue
        qual = fields[5] if len(fields) > 5 else "."
        info = fields[7]
        fmt = fields[8] if len(fields) > 8 else ""
        sample = fields[9] if len(fields) > 9 else ""
        dp = "."
        mq = "."
        bq = qual if qual not in {".", ""} else "."
        info_dict = {}
        for item in info.split(";"):
            if "=" in item:
                k, v = item.split("=", 1)
                info_dict[k] = v
        if "DP" in info_dict:
            dp = info_dict["DP"]
        if "MQ" in info_dict:
            mq = info_dict["MQ"]
        if fmt and sample:
            fmt_dict = dict(zip(fmt.split(":"), sample.split(":")))
            if dp == "." and "DP" in fmt_dict:
                dp = fmt_dict["DP"]
            if mq == "." and "MQ" in fmt_dict:
                mq = fmt_dict["MQ"]
            if bq == "." and "BQ" in fmt_dict:
                bq = fmt_dict["BQ"]
        variant_positions.add((chrom, pos))
        stats_by_pos[(chrom, pos)] = {"BQ": bq, "DP": dp, "MQ": mq}
    return variant_positions, stats_by_pos


def lineage_prediction_from_vcf(vcf_file, snp_db, bed_file):
    lineage_refs = load_lineage_snp_db(snp_db)
    bed_positions = load_bed_positions(bed_file)
    variant_positions, stats_by_pos = parse_vcf_variants(vcf_file)
    if not variant_positions:
        return {"predicted_lineage": "No variant found in MTB reads", "matched_snp_positions": "NA", "representative_bq": "NA", "representative_dp": "NA", "representative_mq": "NA"}
    variant_positions_in_bed = {x for x in variant_positions if x in bed_positions}
    if not variant_positions_in_bed:
        return {"predicted_lineage": "Unpredictable: insufficient SNP evidence", "matched_snp_positions": "NA", "representative_bq": "NA", "representative_dp": "NA", "representative_mq": "NA"}
    results = {}
    matched_positions_per_lineage = {}
    for lineage, ref_snps in lineage_refs.items():
        matched = ref_snps & variant_positions_in_bed
        total = len(ref_snps)
        results[lineage] = round((len(matched) / total) * 100, 2) if total > 0 else 0.0
        matched_positions_per_lineage[lineage] = sorted(matched, key=lambda x: (x[0], x[1]))
    if not results or max(results.values()) == 0:
        return {"predicted_lineage": "Unpredictable: insufficient SNP evidence", "matched_snp_positions": "NA", "representative_bq": "NA", "representative_dp": "NA", "representative_mq": "NA"}
    top_lineage = max(results, key=results.get)
    top_prob = results[top_lineage]
    predicted_set = sorted([lin for lin, pc in results.items() if pc >= 90.0 or pc >= (top_prob - 25.0)])
    predicted = ";".join(predicted_set) if predicted_set else "Unpredictable: insufficient SNP evidence"
    matched_positions = matched_positions_per_lineage.get(top_lineage, [])
    if not matched_positions:
        return {"predicted_lineage": predicted, "matched_snp_positions": "NA", "representative_bq": "NA", "representative_dp": "NA", "representative_mq": "NA"}
    pos_strings, bqs, dps, mqs = [], [], [], []
    for pos in matched_positions:
        pos_strings.append(str(pos[1]))
        st = stats_by_pos.get(pos, {})
        if st.get("BQ", ".") not in [".", "NA", ""]: bqs.append(str(st["BQ"]))
        if st.get("DP", ".") not in [".", "NA", ""]: dps.append(str(st["DP"]))
        if st.get("MQ", ".") not in [".", "NA", ""]: mqs.append(str(st["MQ"]))
    return {"predicted_lineage": predicted, "matched_snp_positions": ",".join(pos_strings) if pos_strings else "NA", "representative_bq": ",".join(bqs) if bqs else "NA", "representative_dp": ",".join(dps) if dps else "NA", "representative_mq": ",".join(mqs) if mqs else "NA"}


def extract_mtb_reads(input_type, original_fastq_files, kraken_out, kraken_report, extract_script, outdir, sample_name):
    ensure_dir(outdir)
    taxids = get_taxids_from_kraken_report(kraken_report)
    tmp_extract_dir = os.path.join(outdir, "tmp_myco_extract")
    ensure_dir(tmp_extract_dir)
    if input_type == "paired":
        myco_r1 = os.path.join(outdir, f"{sample_name}.mycobacterium_reads_R1.fastq")
        myco_r2 = os.path.join(outdir, f"{sample_name}.mycobacterium_reads_R2.fastq")
        open(myco_r1, "w").close(); open(myco_r2, "w").close()
        for taxid in taxids:
            out_r1 = os.path.join(tmp_extract_dir, f"{sample_name}.taxid_{taxid}_R1.fastq")
            out_r2 = os.path.join(tmp_extract_dir, f"{sample_name}.taxid_{taxid}_R2.fastq")
            run_command(["python3", extract_script, "-k", kraken_out, "-r", kraken_report, "-s", original_fastq_files[0], "-s2", original_fastq_files[1], "-t", str(taxid), "--include-children", "--fastq-output", "-o", out_r1, "-o2", out_r2])
            append_file(out_r1, myco_r1); append_file(out_r2, myco_r2)
        dedup_r1 = os.path.join(outdir, f"{sample_name}.mycobacterium_reads_R1.dedup.fastq")
        dedup_r2 = os.path.join(outdir, f"{sample_name}.mycobacterium_reads_R2.dedup.fastq")
        if os.path.getsize(myco_r1) > 0 and os.path.getsize(myco_r2) > 0:
            deduplicate_paired_fastq_by_r1(myco_r1, myco_r2, dedup_r1, dedup_r2)
            os.replace(dedup_r1, myco_r1); os.replace(dedup_r2, myco_r2)
        return [myco_r1, myco_r2]
    myco_se = os.path.join(outdir, f"{sample_name}.mycobacterium_reads.fastq")
    open(myco_se, "w").close()
    for taxid in taxids:
        out_se = os.path.join(tmp_extract_dir, f"{sample_name}.taxid_{taxid}.fastq")
        run_command(["python3", extract_script, "-k", kraken_out, "-r", kraken_report, "-s", original_fastq_files[0], "-t", str(taxid), "--include-children", "--fastq-output", "-o", out_se])
        append_file(out_se, myco_se)
    return [myco_se]


def bwa_to_h37rv(input_type, input_files, ref, outdir, sample_name, threads):
    ensure_dir(outdir)
    sam = os.path.join(outdir, f"{sample_name}.sam")
    bam = os.path.join(outdir, f"{sample_name}.bam")
    sorted_bam = os.path.join(outdir, f"{sample_name}.sorted.bam")
    if input_type == "paired":
        bwa_cmd = f"bwa mem -M -t {threads} -R '@RG\\tID:{sample_name}\\tSM:{sample_name}\\tPL:ILLUMINA' {shlex.quote(ref)} {shlex.quote(input_files[0])} {shlex.quote(input_files[1])} > {shlex.quote(sam)}"
    else:
        bwa_cmd = f"bwa mem -M -t {threads} -R '@RG\\tID:{sample_name}\\tSM:{sample_name}\\tPL:ILLUMINA' {shlex.quote(ref)} {shlex.quote(input_files[0])} > {shlex.quote(sam)}"
    run_command(bwa_cmd, shell=True)
    run_command(["samtools", "view", "-S", "-b", sam, "-o", bam])
    run_command(["samtools", "sort", "-@", str(threads), bam, "-o", sorted_bam])
    run_command(["samtools", "index", sorted_bam])
    return sorted_bam


def bcftools_call(ref, bam, outdir, sample_name):
    ensure_dir(outdir)
    vcf = os.path.join(outdir, f"{sample_name}.vcf")
    vcfgz = f"{vcf}.gz"
    mpileup_cmd = f"bcftools mpileup -Ou -f {shlex.quote(ref)} {shlex.quote(bam)} | bcftools call --ploidy 1 -mv -Ov -o {shlex.quote(vcf)}"
    run_command(mpileup_cmd, shell=True)
    run_command(["bgzip", "-f", vcf])
    run_command(["tabix", "-f", "-p", "vcf", vcfgz])
    return vcfgz


def write_summary_csv(rows, out_csv):
    with open(out_csv, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=SUMMARY_FIELDS)
        writer.writeheader()
        writer.writerows(rows)


def process_sample(entry, args):
    input_type, input_files, sample_name = auto_detect_fastq(entry)
    LOGGER.info("[%s] detected input_type=%s files=%s", sample_name, input_type, ",".join(input_files))
    sample_tmp = os.path.join(args.tmp_dir, sample_name)
    sample_res = os.path.join(args.result_dir, sample_name)
    kraken_dir = os.path.join(sample_tmp, "kraken")
    mtb_dir = os.path.join(sample_res, "Mycobacterial_fastq_files")
    align_dir = os.path.join(sample_tmp, "alignment")
    vcf_dir = os.path.join(sample_tmp, "vcf")
    for d in [sample_tmp, sample_res, kraken_dir, mtb_dir, align_dir, vcf_dir]:
        ensure_dir(d)
    kraken_out = os.path.join(kraken_dir, f"{sample_name}.kraken")
    kraken_report = os.path.join(kraken_dir, f"{sample_name}.report")
    cmd = ["kraken2", "--db", args.kraken_db, "--threads", str(args.threads), "--gzip-compressed", "--report", kraken_report, "--output", kraken_out]
    if args.kraken_memory_mapping:
        cmd.append("--memory-mapping")
    if input_type == "paired":
        cmd.extend(["--paired", input_files[0], input_files[1]])
    else:
        cmd.append(input_files[0])
    run_command(cmd)
    mtb_fastqs = extract_mtb_reads(input_type, input_files, kraken_out, kraken_report, args.extract_script, mtb_dir, sample_name)
    sorted_bam = bwa_to_h37rv(input_type, mtb_fastqs, args.ref, align_dir, sample_name, args.threads)
    vcf_gz = bcftools_call(args.ref, sorted_bam, vcf_dir, sample_name)
    lineage = lineage_prediction_from_vcf(vcf_gz, args.snp_db, args.bed_file)
    return {"Sample_name": sample_name, "SNP_position_Lineage_express_db": lineage["matched_snp_positions"], "Base_quality": lineage["representative_bq"], "Depth": lineage["representative_dp"], "Mapping_Quality": lineage["representative_mq"], "Predicted_lineage": lineage["predicted_lineage"]}


def main():
    p = argparse.ArgumentParser(description="Run lineage prediction from host-removed FASTQ reads (haploid bcftools calling)")
    p.add_argument("--input")
    p.add_argument("--input_list")
    p.add_argument("--input_dir")
    p.add_argument("--ref", default="/home/muthukumarb/drugresistance/data/reference/h37rv.fa")
    p.add_argument("--snp_db", default="/home/muthukumarb/drugresistance/data/database/lineage_snp_updated_au13.tsv")
    p.add_argument("--bed_file", default="/home/muthukumarb/drugresistance/data/database/targeted_modified_regions_au13.bed")
    p.add_argument("--kraken_db", required=True)
    p.add_argument("--extract_script", required=True)
    p.add_argument("--tmp_dir", default="tmp_folder")
    p.add_argument("--result_dir", default="result_folder")
    p.add_argument("--threads", type=int, default=4)
    p.add_argument("--kraken_memory_mapping", action="store_true")
    p.add_argument("--summary_csv", default="final_lineage_summary.csv")
    args = p.parse_args()
    if sum([bool(args.input), bool(args.input_list), bool(args.input_dir)]) != 1:
        raise ValueError("Provide exactly one of --input, --input_list, or --input_dir")
    for label, path in {"ref": args.ref, "snp_db": args.snp_db, "bed_file": args.bed_file, "kraken_db": args.kraken_db, "extract_script": args.extract_script}.items():
        if not os.path.exists(path):
            raise FileNotFoundError(f"{label} not found: {path}")
    ensure_dir(args.tmp_dir); ensure_dir(args.result_dir)
    samples = []
    if args.input_list:
        with open(args.input_list) as fh:
            samples = [x.strip() for x in fh if x.strip()]
    elif args.input:
        samples = [args.input.strip()]
    else:
        seen_entries = set()
        for fname in sorted(os.listdir(args.input_dir)):
            full = os.path.join(args.input_dir, fname)
            if not os.path.isfile(full): continue
            if not any(fname.endswith(ext) for ext in FASTQ_EXTENSIONS): continue
            sample_name = clean_sample_name(fname)
            prefix_path = os.path.join(args.input_dir, sample_name)
            if prefix_path not in seen_entries:
                seen_entries.add(prefix_path)
                samples.append(prefix_path)
    if not samples:
        raise ValueError("No valid samples found to process")
    rows = [process_sample(entry, args) for entry in samples]
    out_csv = args.summary_csv if os.path.isabs(args.summary_csv) else os.path.join(args.result_dir, args.summary_csv)
    write_summary_csv(rows, out_csv)
    LOGGER.info("Wrote summary CSV: %s", out_csv)


if __name__ == "__main__":
    main()
