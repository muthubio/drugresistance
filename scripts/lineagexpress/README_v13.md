# LineageXpress v13 (Optimized)

LineageXpress v13 is an optimized version of v12 with major performance improvements:

## Key improvements over v12

### 1. Targeted resistance variant calling

- Uses a separate `--resistance_bed`
- Avoids full genome variant calling
- Dramatically reduces runtime

### 2. Combined targeted calling

- Merges lineage BED + resistance BED
- Single GATK HaplotypeCaller run
- One VCF used for both lineage and resistance

### 3. Multi-sample parallel execution

- Uses Python multiprocessing
- Multiple samples processed simultaneously

---

## Required inputs

### 1. Lineage BED

```bash
targeted_modified_regions_au13.bed
```

### 2. Resistance BED (NEW)

Example:

```text
NC_000962.3	761155	761565	rpoB
NC_000962.3	2155168	2156340	katG
NC_000962.3	1673425	1673550	inhA_promoter
NC_000962.3	757000	758500	gyrA
NC_000962.3	681000	682500	gyrB
NC_000962.3	4246000	4247600	pncA
```

---

## Example command

```bash
python lineageexpress_new_version13.py \
  --sample_list samples_list.txt \
  --ref_genome /home/seq/genomics/lineage_express/data/reference/h37rv.fa \
  --output_dir results_v13 \
  --snp_file /home/seq/genomics/lineage_express/data/database/lineage_snp_updated_au13.tsv \
  --bed_file /home/seq/genomics/lineage_express/data/database/targeted_modified_regions_au13.bed \
  --resistance_bed /home/seq/genomics/lineage_express/data/database/resistance_targets.bed \
  --mutations_csv /home/seq/genomics/lineage_express/data/database/mutations.csv \
  --threads 8 \
  --max_parallel_samples 4
```

---

## Performance comparison

| Version | Resistance calling | Speed |
|--------|-------------------|------|
| v12 | Full genome | Slow |
| v13 | Targeted BED | Fast ⚡ |

---

## Parallel execution

```bash
--max_parallel_samples 4
```

Runs 4 samples simultaneously.

---

## Output structure

Same as v12.

---

## Recommended workflow

1. Map reads to full H37Rv (no change)
2. Use combined BED for variant calling
3. Perform lineage + resistance from same VCF

---

## Important note

If a SNP is not present in BED:

- It will NOT be called
- It will NOT be used for prediction

So ensure:

- lineage BED covers all lineage SNPs
- resistance BED covers all resistance loci

---

## When to use v13

Use v13 for:

- Large datasets
- HPC execution
- Faster turnaround

---

## Future improvements

- Auto BED generation from SNP database
- Slurm job array support
- Benchmark module
