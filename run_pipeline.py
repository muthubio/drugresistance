#!/usr/bin/env python3
"""
run_pipeline.py
Thin CLI -> pipeline_core.run_pipeline
Now supports:
  --dr-catalog tbdb.dr.json
  --mask-bed tbdb.mask.bed
  --min-alt-af, --min-dp
"""
import argparse
from pipeline_core import (
    PipelineConfig, run_pipeline, INFO, WARNING,
    set_log_level, set_log_file, close_logs
)

def build_parser() -> argparse.ArgumentParser:
    ap = argparse.ArgumentParser(
        description=("Trim Galore -> BWA -> merge -> fix/sort/markdup -> index; "
                     "optional GATK HC + bcftools norm + TB-Profiler-like DR annotation.")
    )
    # IO
    ap.add_argument("--r1", required=True, help="R1 FASTQ(.gz)")
    ap.add_argument("--r2", required=True, help="R2 FASTQ(.gz)")
    ap.add_argument("--ref", required=True, help="Reference FASTA")
    ap.add_argument("--sample", required=True, help="Sample name for @RG:ID and SM")
    ap.add_argument("--out", default=None, help="Output BAM (default: <sample>.bam)")
    ap.add_argument("--tmpdir", default=None, help="Temporary directory (kept if provided)")

    # mapping/trimming
    ap.add_argument("--length", type=int, default=36, help="Minimum read length after trimming (default 36)")
    ap.add_argument("--threads", type=int, default=1, help="Threads for BWA/Samtools/GATK")
    ap.add_argument("--mem", default="768M", help="Per-sort memory for samtools (default 768M)")
    ap.add_argument("--fastqc", action="store_true", help="Run FastQC via Trim Galore")
    ap.add_argument("--skip-singles", action="store_true", help="Skip mapping unpaired reads")

    # calling
    ap.add_argument("--call-variants", action="store_true", help="Run GATK HaplotypeCaller after BAM creation")
    ap.add_argument("--emit-gvcf", action="store_true", help="Emit gVCF (HaplotypeCaller -ERC GVCF)")
    ap.add_argument("--intervals", default=None, help="Intervals list/BED for HaplotypeCaller (-L)")
    ap.add_argument("--normalize", action="store_true",
                    help="Normalize VCF with bcftools norm (-f REF, split multiallelics)")
    ap.add_argument("--gatk-extra", default=None, help="Extra args to append to HaplotypeCaller (quoted string)")

    # DR annotation
    ap.add_argument("--dr-catalog", default=None, help="Drug-resistance JSON (e.g., tbdb.dr.json)")
    ap.add_argument("--mask-bed", default=None, help="BED of regions to exclude before DR annotation")
    ap.add_argument("--min-alt-af", type=float, default=0.0, help="Minimum ALT allele fraction to count")
    ap.add_argument("--min-dp", type=int, default=0, help="Minimum depth to count")

    # logging
    ap.add_argument("--logging", choices=["DEBUG", "INFO", "WARNING", "ERROR"],
                    default="INFO", help="Console log level (default: INFO)")
    ap.add_argument("--log-file", dest="log_file", default=None,
                    help="Also append logs to this file")
    return ap

def main():
    ap = build_parser()
    args = ap.parse_args()

    set_log_level(args.logging)
    set_log_file(args.log_file)

    if args.normalize and not args.call_variants:
        INFO("--normalize has no effect without --call-variants; ignoring.")
    if (args.dr_catalog or args.mask_bed) and not args.call_variants:
        WARNING("--dr-catalog / --mask-bed provided but --call-variants is off; no VCF => DR/masking skipped.")

    cfg = PipelineConfig(
        r1=args.r1, r2=args.r2, ref=args.ref, sample=args.sample,
        out_bam=args.out, tmpdir=args.tmpdir,
        length=args.length, threads=args.threads, mem=args.mem,
        fastqc=args.fastqc, skip_singles=args.skip_singles,
        call_variants=args.call_variants, emit_gvcf=args.emit_gvcf,
        intervals=args.intervals, normalize=args.normalize, gatk_extra=args.gatk_extra,
        dr_catalog=args.dr_catalog, min_alt_af=args.min_alt_af, min_dp=args.min_dp,
        mask_bed=args.mask_bed
    )

    try:
        outputs = run_pipeline(cfg)
        INFO("=== Outputs ===")
        for k, v in outputs.items():
            INFO(f"{k:12s}: {v}")
    finally:
        close_logs()

if __name__ == "__main__":
    main()
