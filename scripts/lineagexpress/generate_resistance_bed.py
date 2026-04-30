#!/usr/bin/env python3
"""
generate_resistance_bed.py

Generate a resistance target BED file for LineageXpress v13 from:
  1. a mutations.csv file, usually WHO/TB-Profiler style
  2. a gene coordinate table for H37Rv

This is safer than manually typing resistance BED coordinates.

Why needed?
----------
LineageXpress v13 uses targeted resistance variant calling:

  gatk HaplotypeCaller -L resistance_targets.bed

If a resistance gene/promoter is missing from the BED file, variants in that region
will not be called. This script helps build the BED systematically.

Expected mutations.csv columns
------------------------------
The script auto-detects common columns:
  Gene / gene / locus / Rv_locus

Expected gene coordinate TSV columns
------------------------------------
Required columns:
  gene    chrom    start    end

Optional columns:
  strand
  promoter_upstream
  notes

Coordinate convention
---------------------
The coordinate TSV may use 1-based inclusive gene coordinates.
The output BED is 0-based start, half-open end, as required by BED.

Example gene coordinate TSV:
----------------------------
gene	chrom	start	end	strand	promoter_upstream
rpoB	NC_000962.3	761155	763325	+	0
katG	NC_000962.3	2153889	2156111	-	0
inhA	NC_000962.3	1674202	1675011	+	200
fabG1	NC_000962.3	1673425	1674201	+	200
gyrA	NC_000962.3	7302	9818	+	0
gyrB	NC_000962.3	5240	7267	+	0
pncA	NC_000962.3	2288681	2289241	+	200
embB	NC_000962.3	4246514	4249810	+	0
rrs	NC_000962.3	1471846	1473382	+	0
eis	NC_000962.3	2714124	2714817	-	200

Usage
-----
python generate_resistance_bed.py \
  --mutations_csv /home/seq/genomics/lineage_express/data/database/mutations.csv \
  --gene_coordinates resistance_gene_coordinates_h37rv.tsv \
  --out resistance_targets.bed \
  --padding 50

Recommended output path for v13:
--------------------------------
/home/seq/genomics/lineage_express/data/database/resistance_targets.bed
"""

import argparse
import csv
from pathlib import Path
from typing import Dict, List, Set, Tuple


def norm(s):
    return "" if s is None else str(s).strip()


def norm_col(s):
    return norm(s).lower().replace(" ", "_").replace("-", "_")


def detect_delimiter(path: str) -> str:
    with open(path, "r", encoding="utf-8", errors="ignore") as f:
        first = f.readline()
    return "\t" if "\t" in first else ","


def read_table(path: str) -> List[Dict[str, str]]:
    delim = detect_delimiter(path)
    with open(path, "r", encoding="utf-8", errors="ignore") as f:
        reader = csv.DictReader(f, delimiter=delim)
        if reader.fieldnames is None:
            raise SystemExit(f"No header found in {path}")
        reader.fieldnames = [norm_col(x) for x in reader.fieldnames]
        rows = []
        for row in reader:
            rows.append({norm_col(k): norm(v) for k, v in row.items()})
    return rows


def pick_gene_column(rows: List[Dict[str, str]]) -> str:
    candidates = [
        "gene", "genes", "locus", "rv_locus", "rvlocus", "gene_name",
        "target_gene", "gene_or_locus", "gene_id"
    ]
    if not rows:
        raise SystemExit("mutations_csv has no rows")
    cols = set(rows[0].keys())
    for c in candidates:
        if c in cols:
            return c
    raise SystemExit(
        "Could not detect gene column in mutations_csv. Expected one of: "
        + ", ".join(candidates)
    )


def split_gene_field(value: str) -> List[str]:
    value = norm(value)
    if not value:
        return []
    for sep in [";", "|", "/", ","]:
        value = value.replace(sep, ";")
    genes = []
    for g in value.split(";"):
        g = g.strip()
        if g:
            genes.append(g)
    return genes


def load_mutation_genes(mutations_csv: str) -> Set[str]:
    rows = read_table(mutations_csv)
    gene_col = pick_gene_column(rows)
    genes = set()
    for row in rows:
        for g in split_gene_field(row.get(gene_col, "")):
            genes.add(g)
    return genes


def load_gene_coordinates(path: str) -> Dict[str, Dict[str, str]]:
    rows = read_table(path)
    required = {"gene", "chrom", "start", "end"}
    cols = set(rows[0].keys()) if rows else set()
    missing = required - cols
    if missing:
        raise SystemExit(f"gene coordinate file missing columns: {', '.join(sorted(missing))}")

    coords = {}
    for row in rows:
        gene = row["gene"]
        if not gene:
            continue
        coords[gene.lower()] = row
    return coords


def make_bed_interval(row: Dict[str, str], padding: int) -> Tuple[str, int, int, str]:
    gene = row["gene"]
    chrom = row["chrom"]
    try:
        start_1based = int(float(row["start"]))
        end_1based = int(float(row["end"]))
    except ValueError:
        raise SystemExit(f"Invalid coordinates for gene {gene}: start={row['start']} end={row['end']}")

    strand = row.get("strand", "")
    promoter_upstream = row.get("promoter_upstream", "0") or "0"
    try:
        promoter = int(float(promoter_upstream))
    except ValueError:
        promoter = 0

    start = min(start_1based, end_1based)
    end = max(start_1based, end_1based)

    # Convert to BED: 0-based start, half-open end.
    bed_start = max(0, start - 1 - padding)
    bed_end = end + padding

    # Add promoter region directionally if requested.
    if promoter > 0:
        if strand == "-":
            bed_end = end + promoter + padding
        else:
            bed_start = max(0, start - 1 - promoter - padding)

    return chrom, bed_start, bed_end, gene


def main():
    p = argparse.ArgumentParser(
        description="Generate resistance_targets.bed from mutations.csv and H37Rv gene coordinate table"
    )
    p.add_argument("--mutations_csv", required=True, help="WHO/TB-Profiler style mutations.csv")
    p.add_argument("--gene_coordinates", required=True, help="TSV/CSV with gene,chrom,start,end,strand,promoter_upstream")
    p.add_argument("--out", required=True, help="Output BED file")
    p.add_argument("--padding", type=int, default=50, help="Extra bases added on both sides of target interval")
    p.add_argument("--include_all_coordinates", action="store_true", help="Ignore mutations.csv genes and output all genes in coordinate file")
    args = p.parse_args()

    if not Path(args.mutations_csv).exists():
        raise SystemExit(f"mutations_csv not found: {args.mutations_csv}")
    if not Path(args.gene_coordinates).exists():
        raise SystemExit(f"gene_coordinates not found: {args.gene_coordinates}")

    coords = load_gene_coordinates(args.gene_coordinates)
    mutation_genes = load_mutation_genes(args.mutations_csv)

    if args.include_all_coordinates:
        selected_keys = sorted(coords.keys())
    else:
        selected_keys = []
        missing = []
        for gene in sorted(mutation_genes):
            key = gene.lower()
            if key in coords:
                selected_keys.append(key)
            else:
                missing.append(gene)

        if missing:
            print("[WARNING] Genes present in mutations.csv but missing in coordinate table:")
            for g in missing:
                print(f"  - {g}")

    intervals = []
    for key in selected_keys:
        intervals.append(make_bed_interval(coords[key], args.padding))

    # Sort and merge exact overlapping intervals per chromosome.
    intervals.sort(key=lambda x: (x[0], x[1], x[2], x[3]))

    out_path = Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, "w") as out:
        for chrom, start, end, gene in intervals:
            out.write(f"{chrom}\t{start}\t{end}\t{gene}\n")

    print(f"[DONE] Wrote BED: {out_path}")
    print(f"[INFO] Genes in mutations.csv: {len(mutation_genes)}")
    print(f"[INFO] Genes written to BED: {len(intervals)}")
    print(f"[INFO] Padding: {args.padding} bp")


if __name__ == "__main__":
    main()
