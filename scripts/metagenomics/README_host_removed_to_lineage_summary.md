# Host-Removed FASTQ to Lineage Summary Pipeline

## Overview
Pipeline for MTBC lineage prediction starting from host-removed metagenomics reads:

Host-removed FASTQ -> Kraken2 -> MTBC read extraction -> BWA mapping -> BAM -> bcftools variants -> lineage prediction -> final summary CSV

---

## Requirements

- bwa
- samtools
- bcftools
- kraken2
- KrakenTools (`extract_kraken_reads.py`)
- bgzip
- tabix

---

## Usage

### Single sample

```bash
python host_removed_to_lineage_summary.py \
  --input SRR12569816_host_removed_R1.fastq.gz,SRR12569816_host_removed_R2.fastq.gz \
  --kraken_db /path/to/kraken_db \
  --extract_script /path/to/extract_kraken_reads.py \
  --threads 4
```

### Multiple samples from a list

```bash
python host_removed_to_lineage_summary.py \
  --input_list samples.txt \
  --kraken_db /path/to/kraken_db \
  --extract_script /path/to/extract_kraken_reads.py \
  --threads 4
```

```bash
 python host_removed_to_lineage_summary.py \
  --input_list sample_list.txt \
  --kraken_db /home/muthukumarb/metagenomics_tb/data/kraken_gtdb \
  --extract_script /home/muthukumarb/krakendb/KrakenTools/extract_kraken_reads.py \
  --kraken_memory_mapping \
  --threads 12 \
  --result_dir /home/muthukumarb/metagenomics_tb/results/mtb_lineage_express_testing/test8 \
  --summary_csv /home/muthukumarb/metagenomics_tb/results/mtb_lineage_express_testing/test8.csv
```

### Automatic detection from a directory

```bash
python host_removed_to_lineage_summary.py \
  --input_dir /path/to/host_removed_fastq \
  --kraken_db /path/to/kraken_db \
  --extract_script /path/to/extract_kraken_reads.py \
  --threads 4
```

---

## Input formats supported

- Paired host-removed FASTQ
- Single-end host-removed FASTQ
- Sample prefix path (auto-detect)
- Input list file for multiple samples

---

## Output

- Mycobacterial extracted FASTQ
- Sorted BAM
- VCF (`.vcf.gz` + index)
- Final summary CSV

Final summary columns:

- `Sample_name`
- `SNP_position_Lineage_express_db`
- `Base_quality`
- `Depth`
- `Mapping_Quality`
- `Predicted_lineage`

---

## Notes

- Designed for MTBC lineage prediction from **host-removed reads**
- Starts after host removal, so FastQC / trimming / human filtering are not repeated
- Uses Kraken2 + taxonomy-based MTB read extraction before alignment to H37Rv
- Uses haploid variant calling
- Supports single or multiple samples in one run
