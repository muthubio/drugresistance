# Metagenomics BWA Mapping Pipeline

## Overview
Pipeline for MTBC lineage prediction from metagenomics reads:

FASTQ -> BWA mapping -> BAM -> bcftools variants -> lineage prediction

## Requirements
- bwa
- samtools
- bcftools

## Usage
```bash
python Metagenomics_bwa_mapping_script.py \
  --sample_list samples.txt \
  --ref_genome data/h37rv.fa \
  --snp_file data/lineage_snp_updated_au13.tsv \
  --threads 4
```

## Input formats supported
- Paired FASTQ
- Single FASTQ
- BAM
- VCF

## Output
- Sorted BAM
- VCF
- Lineage result file

## Notes
- Designed for MTBC lineage prediction
- Uses haploid variant calling
