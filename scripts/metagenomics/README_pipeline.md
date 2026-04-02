# Metagenomics BWA Mapping Pipeline

## Overview
Pipeline for MTBC lineage prediction from metagenomics reads:

FASTQ -> FastQC -> Trim Galore -> Host Removal -> Kraken2 -> BWA mapping -> BAM -> bcftools variants -> lineage prediction

---

## Requirements

- bwa
- samtools
- bcftools
- fastqc
- trim_galore
- kraken2
- KrakenTools (extract_kraken_reads.py)

---

## Usage

```bash
python mtb_meta_lineage_monitor_ver4.py \
  --input_list samples.txt \
  --human_ref GRCh38.fa \
  --extract_script extract_kraken_reads.py \
  --threads 4
```

---

## Input formats supported

- Paired FASTQ
- Single FASTQ
- Sample prefix (auto-detect)

---

## Output

- Host removed FASTQ
- Mycobacterial extracted FASTQ
- Sorted BAM
- VCF (bgzip + tabix indexed)
- Lineage result file
- Final summary TSV
- Per-step logs
- JSON monitoring state file

---

## Notes

- Designed for MTBC lineage prediction from metagenomics data
- Uses haploid variant calling
- Supports monitoring and failure tracking
- Automatically detects paired/single-end reads
- Uses Kraken2 + taxonomy-based read extraction before lineage prediction
