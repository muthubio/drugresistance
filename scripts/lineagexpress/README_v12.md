# LineageXpress v12

LineageXpress v12 predicts MTBC lineage and optionally performs drug-resistance mutation matching using a WHO/TB-Profiler-style mutation catalogue.

## Main features

- Supports FASTQ, BAM, and VCF input through `--sample_list`.
- Maps FASTQ reads to the full H37Rv reference using BWA MEM.
- Calls lineage variants using GATK HaplotypeCaller restricted to `--bed_file`.
- Performs resistance prediction when `--mutations_csv` is provided.
- Generates lineage result text files and an infographic PDF report.
- Can optionally run Delly for structural variant calling unless `--skip_delly` is used.

## Important note about v12 resistance calling

In v12, lineage calling is targeted using `--bed_file`, but resistance GATK HaplotypeCaller is run on the full genome. This can be slow for many FASTQ/BAM samples.

For faster production runs, use v13, which supports a separate `--resistance_bed` and targeted resistance variant calling.

## Example command used on server

```bash
python lineageexpress_new_version12.py \
  --sample_list samples_list.txt \
  --ref_genome /home/seq/genomics/lineage_express/data/reference/h37rv.fa \
  --output_dir retest \
  --snp_file /home/seq/genomics/lineage_express/data/database/lineage_snp_updated_au13.tsv \
  --bed_file /home/seq/genomics/lineage_express/data/database/targeted_modified_regions_au13.bed \
  --mutations_csv /home/seq/genomics/lineage_express/data/database/mutations.csv
```

## Recommended v12 command with more threads

```bash
python lineageexpress_new_version12.py \
  --sample_list samples_list.txt \
  --ref_genome /home/seq/genomics/lineage_express/data/reference/h37rv.fa \
  --output_dir retest \
  --snp_file /home/seq/genomics/lineage_express/data/database/lineage_snp_updated_au13.tsv \
  --bed_file /home/seq/genomics/lineage_express/data/database/targeted_modified_regions_au13.bed \
  --mutations_csv /home/seq/genomics/lineage_express/data/database/mutations.csv \
  --threads 8
```

## Sample list format

One sample per line. Accepted formats:

```text
/path/sample_R1.fastq.gz,/path/sample_R2.fastq.gz
/path/sample.fastq.gz
/path/sample.sorted.bam
/path/sample.vcf
sample_id
```

If only a bare sample ID is provided, the script searches current directory, `sample_input/`, and `sample_data/`.

## Main outputs

For each sample, outputs are written under:

```text
<output_dir>/<sample_id>/
```

Common output files:

```text
<sample_id>_lineage_result.txt
matched_mutations.csv
mutations_summary.csv
<sample_id>_final_prediction.txt
<sample_id>_summary.csv
<sample_id>_infographic_report.pdf
```

## Dependencies

- Python 3
- pandas
- reportlab
- bwa
- samtools
- gatk
- bcftools
- snpEff
- delly, optional if using SV calling

## When to use v12

Use v12 when you want the older stable workflow and do not mind full-genome resistance calling.

Use v13 when you want faster multi-sample execution and targeted resistance calling.
