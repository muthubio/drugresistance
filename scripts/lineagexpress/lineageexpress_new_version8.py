#!/usr/bin/env python3
import os
import re
import shlex
import subprocess
import argparse
import logging
import time
from collections import defaultdict
from pathlib import Path
import pandas as pd

from reportlab.lib import colors
from reportlab.lib.pagesizes import A3, landscape
from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
from reportlab.lib.enums import TA_LEFT, TA_CENTER
from reportlab.lib.units import mm
from reportlab.platypus import (
    SimpleDocTemplate, Paragraph, Spacer, Table, TableStyle,
    Flowable, KeepTogether
)

logging.basicConfig(format='%(asctime)s - %(levelname)s - %(message)s', level=logging.INFO)

BASE_DIR = Path(__file__).resolve().parent
DATA_DIR = BASE_DIR / "data"
DEFAULT_REF       = str(DATA_DIR / "h37rv.fa")
DEFAULT_SNP_FILE  = str(DATA_DIR / "lineage_snp_updated_au13.tsv")
DEFAULT_BED_FILE  = str(DATA_DIR / "targeted_modified_regions_au13.bed")
DEFAULT_SNPEFF_DB = "Mycobacterium_tuberculosis_h37rv"

SAMPLE_INPUT_DIR = BASE_DIR / "sample_input"
SAMPLE_DATA_DIR  = BASE_DIR / "sample_data"

def run_cmd(cmd: str) -> bool:
    logging.info(f"Running: {cmd}")
    try:
        subprocess.run(cmd, check=True, shell=True)
        return True
    except subprocess.CalledProcessError as e:
        logging.error(f"Command failed: {cmd}\n{e}")
        return False

def q(x) -> str:
    return shlex.quote(str(x))

def _is_fastq_name(name: str) -> bool:
    low = name.lower()
    return low.endswith((".fastq.gz", ".fq.gz", ".fastq", ".fq"))

def _is_bam_name(name: str) -> bool:
    return name.lower().endswith(".bam")

def _is_vcf_name(name: str) -> bool:
    low = name.lower()
    return low.endswith((".vcf", ".vcf.gz"))

def _classify_fastq_files(paths):
    fq_paths = sorted([Path(p) for p in paths if _is_fastq_name(Path(p).name)], key=lambda x: x.name)
    if not fq_paths:
        return "unknown", []
    by_name = {p.name: p for p in fq_paths}
    explicit_pairs = [
        ("_R1.fastq.gz", "_R2.fastq.gz"), ("_1.fastq.gz", "_2.fastq.gz"), (".1.fastq.gz", ".2.fastq.gz"),
        ("_R1.fq.gz", "_R2.fq.gz"), ("_1.fq.gz", "_2.fq.gz"), (".1.fq.gz", ".2.fq.gz"),
        ("_R1.fastq", "_R2.fastq"), ("_1.fastq", "_2.fastq"), (".1.fastq", ".2.fastq"),
        ("_R1.fq", "_R2.fq"), ("_1.fq", "_2.fq"), (".1.fq", ".2.fq"),
    ]
    for p in fq_paths:
        nm = p.name
        for s1, s2 in explicit_pairs:
            if nm.endswith(s1):
                mate = nm[:-len(s1)] + s2
                if mate in by_name:
                    return "paired", [str(p), str(by_name[mate])]
    if len(fq_paths) == 1:
        return "single", [str(fq_paths[0])]
    if len(fq_paths) == 2:
        return "paired", [str(fq_paths[0]), str(fq_paths[1])]
    return "unknown", []

def _find_in_dir(dirpath: Path, sid: str):
    pats = [
        (f"{sid}_1.fastq.gz",  f"{sid}_2.fastq.gz"),
        (f"{sid}_R1.fastq.gz", f"{sid}_R2.fastq.gz"),
        (f"{sid}.1.fastq.gz",  f"{sid}.2.fastq.gz"),
        (f"{sid}_1.fastq",     f"{sid}_2.fastq"),
        (f"{sid}_R1.fastq",    f"{sid}_R2.fastq"),
        (f"{sid}.1.fastq",     f"{sid}.2.fastq"),
        (f"{sid}_1.fq.gz",     f"{sid}_2.fq.gz"),
        (f"{sid}_R1.fq.gz",    f"{sid}_R2.fq.gz"),
        (f"{sid}.1.fq.gz",     f"{sid}.2.fq.gz"),
        (f"{sid}_1.fq",        f"{sid}_2.fq"),
        (f"{sid}_R1.fq",       f"{sid}_R2.fq"),
        (f"{sid}.1.fq",        f"{sid}.2.fq"),
    ]
    for r1, r2 in pats:
        r1p, r2p = dirpath / r1, dirpath / r2
        if r1p.exists() and r2p.exists():
            return "paired", [str(r1p), str(r2p)]
    for ext in (".fastq.gz", ".fq.gz", ".fastq", ".fq"):
        se = dirpath / f"{sid}{ext}"
        if se.exists():
            return "single", [str(se)]
    try:
        entries = [p for p in dirpath.iterdir() if p.is_file()]
    except FileNotFoundError:
        entries = []
    bam_files = [str(p) for p in entries if _is_bam_name(p.name)]
    if len(bam_files) == 1:
        return "bam", bam_files
    vcf_files = [str(p) for p in entries if _is_vcf_name(p.name)]
    if len(vcf_files) == 1:
        return "vcf", vcf_files
    kind, files = _classify_fastq_files(entries)
    if kind != "unknown":
        return kind, files
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
        return _find_in_dir(p, p.name)
    sid = s
    for base in [Path.cwd(), SAMPLE_INPUT_DIR, SAMPLE_DATA_DIR]:
        kind, files = _find_in_dir(base, sid)
        if kind != "unknown":
            return kind, files
    return "unknown", []

def load_lineage_snp_file(snp_file):
    lineage_snps = defaultdict(set)
    with open(snp_file, 'r') as f:
        for line in f:
            if line.strip() and not line.startswith('#'):
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
    return lineage_snps

def pct(matched, total):
    return round((matched / total) * 100, 2) if total > 0 else 0.0

CONTIG_ALIASES = {"Chromosome", "NC_000962.3"}

def contig_alias_to_nc(chrom: str) -> str:
    return "NC_000962.3" if chrom in CONTIG_ALIASES else chrom

def rewrite_bed_contig(bed_in: str, bed_out: str, target: str = "NC_000962.3"):
    with open(bed_in) as fi, open(bed_out, "w") as fo:
        for ln in fi:
            if not ln.strip() or ln.startswith("#"):
                fo.write(ln)
                continue
            parts = ln.rstrip("\n").split("\t")
            if parts:
                parts[0] = target if parts[0] in CONTIG_ALIASES else parts[0]
            fo.write("\t".join(parts) + "\n")

def lineage_snps_to_target(lineage_snps: dict, target: str = "NC_000962.3"):
    converted = {}
    for lin, sset in lineage_snps.items():
        converted[lin] = {(target if c in CONTIG_ALIASES else c, pos) for (c, pos) in sset}
    return converted

def write_chrmap(tmp_map_path: str, to_chr: bool = True):
    with open(tmp_map_path, "w") as fh:
        if to_chr:
            fh.write("NC_000962.3\tChromosome\n")
        else:
            fh.write("Chromosome\tNC_000962.3\n")

AA3_TO_1 = {
    "Ala":"A","Arg":"R","Asn":"N","Asp":"D","Cys":"C","Gln":"Q","Glu":"E","Gly":"G","His":"H",
    "Ile":"I","Leu":"L","Lys":"K","Met":"M","Phe":"F","Pro":"P","Ser":"S","Thr":"T","Trp":"W",
    "Tyr":"Y","Val":"V","Ter":"*","Stop":"*"
}

def _n(s):
    return "" if pd.isna(s) else str(s).strip()

def _g(s):
    return _n(s).lower()

def _rv(s):
    return _n(s).lower()

def _c(s):
    return _n(s).replace(" ", "")

def to_short_protein(hgvs_p: str) -> str:
    s = _n(hgvs_p)
    if not s:
        return ""
    if s.startswith("p."):
        m = re.fullmatch(r"p\.([A-Z][a-z]{2})(\d+)([A-Z][a-z]{2}|\*)", s)
        if m:
            a1, pos, a2 = m.groups()
            a1 = AA3_TO_1.get(a1, "")
            a2 = "*" if a2 in ("*", "Ter", "Stop") else AA3_TO_1.get(a2, "")
            if a1 and a2:
                return f"{a1}{pos}{a2}"
    m = re.search(r"([A-Z\*]\d+[A-Z\*])", s)
    return m.group(1) if m else ""

def parse_snpeff_vcf_to_df(ann_vcf_path: str) -> pd.DataFrame:
    rows = []
    with open(ann_vcf_path, "r", encoding="utf-8", errors="ignore") as fh:
        for line in fh:
            if not line or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t", 8)
            if len(parts) < 8:
                continue
            chrom, pos, _id, ref, alt, _qual, _filter, info = parts[:8]
            ann_val = None
            for kv in info.split(";"):
                if kv.startswith("ANN="):
                    ann_val = kv[4:]
                    break
            if not ann_val:
                continue
            first_ann = ann_val.split(",", 1)[0]
            cols = first_ann.split("|")
            if len(cols) < 16:
                cols += [""] * (16 - len(cols))
            rows.append([
                int(pos), ref, alt,
                cols[1], cols[2], cols[3], cols[4], cols[9], cols[10],
            ])
    return pd.DataFrame(rows, columns=[
        "Position", "REF", "ALT", "Annotation", "Impact", "Gene", "Rv_locus", "NA_change", "AA_change"
    ])

def conf_is_resistant(conf: str) -> bool:
    if not isinstance(conf, str):
        return False
    s = conf.strip().lower()
    if "not assoc" in s or "not associated" in s:
        return False
    return ("assoc" in s) or ("associated" in s) or ("confers" in s) or ("causes" in s) or ("likely" in s)

def match_catalogue(ann_csv: str, mutations_csv: str, out_all: str, out_res: str):
    ann = pd.read_csv(ann_csv)
    ann["Gene_norm"] = ann.get("Gene", "").map(_g)
    ann["Rv_norm"]   = ann.get("Rv_locus", pd.Series([""] * len(ann))).map(_rv)
    ann["cDNA_norm"] = ann.get("NA_change", pd.Series([""] * len(ann))).map(_c)
    ann["AA_short"]  = ann.get("AA_change", pd.Series([""] * len(ann))).map(to_short_protein)
    ann = ann[(ann["Gene_norm"] != "") | (ann["Rv_norm"] != "")].copy()

    muts = pd.read_csv(mutations_csv)
    if "Gene" not in muts.columns:
        muts["Gene"] = ""
    if "Mutation" not in muts.columns:
        muts["Mutation"] = ""
    if "original_mutation" not in muts.columns:
        muts["original_mutation"] = ""
    muts["Gene_norm"] = muts["Gene"].map(_g)
    muts["mut_pick"] = muts["Mutation"].fillna(muts["original_mutation"])
    muts["cDNA_norm"] = muts["mut_pick"].map(_c)
    muts["AA_short"] = muts["mut_pick"].map(to_short_protein)
    muts["Gene_or_rv_norm"] = muts["Gene"].map(_rv)

    matches = []
    if not ann.empty and not muts.empty:
        t1 = ann.merge(muts, how="inner", left_on=["Gene_norm", "cDNA_norm"], right_on=["Gene_norm", "cDNA_norm"], suffixes=("_ann", "_cat"))
        if not t1.empty:
            t1["_tier"] = "T1_gene+cDNA"
            matches.append(t1)
    if not ann.empty and ann["Rv_norm"].ne("").any():
        t2 = ann.merge(muts, how="inner", left_on=["Rv_norm", "cDNA_norm"], right_on=["Gene_or_rv_norm", "cDNA_norm"], suffixes=("_ann", "_cat"))
        if not t2.empty:
            t2["_tier"] = "T2_rv+cDNA"
            matches.append(t2)

    ann_p = ann[ann["AA_short"] != ""].copy()
    muts_p = muts[muts["AA_short"] != ""].copy()
    if not ann_p.empty and not muts_p.empty:
        t3 = ann_p.merge(muts_p, how="inner", left_on=["Gene_norm", "AA_short"], right_on=["Gene_norm", "AA_short"], suffixes=("_ann", "_cat"))
        if not t3.empty:
            t3["_tier"] = "T3_gene+AA"
            matches.append(t3)
        if ann["Rv_norm"].ne("").any():
            t3b = ann_p.merge(muts_p, how="inner", left_on=["Rv_norm", "AA_short"], right_on=["Gene_or_rv_norm", "AA_short"], suffixes=("_ann", "_cat"))
            if not t3b.empty:
                t3b["_tier"] = "T3b_rv+AA"
                matches.append(t3b)

    if not matches:
        pd.DataFrame().to_csv(out_all, index=False)
        pd.DataFrame().to_csv(out_res, index=False)
        return pd.DataFrame(), pd.DataFrame()

    merged = pd.concat(matches, ignore_index=True, sort=False)
    if "Gene" not in merged.columns:
        if "Gene_ann" in merged.columns:
            merged["Gene"] = merged["Gene_ann"]
        elif "Gene_cat" in merged.columns:
            merged["Gene"] = merged["Gene_cat"]
        else:
            merged["Gene"] = ""

    keep_context = [c for c in ["Position", "REF", "ALT", "Annotation", "Impact", "Gene", "Rv_locus", "NA_change", "AA_change", "_tier"] if c in merged.columns]
    keep_catalog = [c for c in ["Gene", "Mutation", "original_mutation", "type", "drug", "confidence", "source", "comment"] if c in merged.columns]
    cols_all = list(dict.fromkeys(keep_context + keep_catalog))
    out_all_df = merged[cols_all].copy()

    mutkey = out_all_df.get("Mutation").fillna(out_all_df.get("original_mutation"))
    dedup_key = pd.DataFrame({
        "Gene": out_all_df.get("Gene"),
        "MutKey": mutkey,
        "drug": out_all_df.get("drug"),
        "tier": out_all_df.get("_tier"),
        "Position": out_all_df.get("Position"),
    })
    out_all_df = out_all_df.loc[~dedup_key.duplicated()].reset_index(drop=True)

    if "Position" in out_all_df.columns:
        try:
            out_all_df["Position"] = pd.to_numeric(out_all_df["Position"], errors="coerce")
        except Exception:
            pass

    sort_order = [c for c in ["Gene", "Position", "drug", "_tier"] if c in out_all_df.columns]
    if sort_order:
        out_all_df = out_all_df.sort_values(sort_order, kind="mergesort")
    out_all_df.to_csv(out_all, index=False)
    out_res_df = out_all_df.copy()
    out_res_df.to_csv(out_res, index=False)
    return out_all_df, out_res_df

def write_drug_summary(matched_res_path: str, out_csv: str):
    try:
        df = pd.read_csv(matched_res_path)
    except Exception:
        return
    if df.empty:
        sensitive_summary = pd.DataFrame({
            "drug": ["pan-susceptible"],
            "Gene": ["N/A"],
            "Mutation": ["No resistance mutations"],
            "confidence": ["Sensitive"],
            "n": [1]
        })
        sensitive_summary.to_csv(out_csv, index=False)
        return
    for c in ["drug", "Gene", "Mutation", "original_mutation", "confidence"]:
        if c not in df.columns:
            df[c] = ""
    df["MutKey"] = df["Mutation"].fillna(df["original_mutation"])
    summary = (
        df.groupby(["drug", "Gene", "MutKey", "confidence"], dropna=False)
        .size().reset_index(name="n")
        .sort_values(["drug", "Gene", "MutKey"])
    )
    summary = summary.rename(columns={"MutKey": "Mutation"})
    summary.to_csv(out_csv, index=False)

def _norm(s):
    return "" if pd.isna(s) else str(s).strip().lower()

RIF_NAMES = {"rifampicin", "rifampin"}
INH_NAMES = {"isoniazid", "inh"}
FQ_NAMES  = {"levofloxacin", "moxifloxacin", "ofloxacin", "gatifloxacin"}
GA_EXTRA  = {"bedaquiline", "linezolid"}

def _present(df, name_set):
    if "confidence" not in df.columns:
        return False
    resistant_df = df[df["confidence"].map(conf_is_resistant)]
    return resistant_df["drug"].map(_norm).isin(name_set).any()

def _evidence_rows(df, name_set, cols=("drug", "Gene", "Mutation", "original_mutation", "confidence")):
    if "confidence" not in df.columns:
        return pd.DataFrame()
    resistant_df = df[df["confidence"].map(conf_is_resistant)].copy()
    cols = [c for c in cols if c in resistant_df.columns]
    return (
        resistant_df[resistant_df["drug"].map(_norm).isin(name_set)][cols]
        .fillna("").drop_duplicates().reset_index(drop=True)
    )

def predict_final_category(matched_res_csv: str, out_txt: str, out_csv: str) -> str:
    if not Path(matched_res_csv).exists():
        return "No data"
    df = pd.read_csv(matched_res_csv)
    if df.empty:
        cat = "Sensitive Strain"
        Path(out_txt).write_text(cat + "\n", encoding="utf-8")
        pd.DataFrame({"category": [cat]}).to_csv(out_csv, index=False)
        return cat
    if "drug" not in df.columns:
        raise SystemExit("matched_mutations.csv has no 'drug' column")

    has_rif = _present(df, RIF_NAMES)
    has_inh = _present(df, INH_NAMES)
    has_fq  = _present(df, FQ_NAMES)
    has_ga  = _present(df, GA_EXTRA)

    if (has_rif or has_inh) and has_fq and has_ga:
        category = "XDR-TB"
    elif (has_rif and has_inh) and has_fq:
        category = "pre-XDR-TB"
    elif has_rif and has_fq:
        category = "pre-XDR-TB"
    elif has_rif and has_inh:
        category = "MDR-TB"
    elif has_rif:
        category = "RR-TB"
    elif has_inh and not has_rif:
        category = "Hr-TB"
    else:
        category = "Sensitive Strain"

    Path(out_txt).write_text(f"Final category: {category}\n", encoding="utf-8")
    pd.DataFrame([{
        "RR_any_rif": bool(has_rif),
        "Hr_any_inh": bool(has_inh),
        "FQ_any": bool(has_fq),
        "GroupA_additional_any": bool(has_ga),
        "category": category
    }]).to_csv(out_csv, index=False)
    return category

class PillBox(Flowable):
    def __init__(self, label, bg, width=44*mm, height=14*mm, font_size=9, text_color=colors.white):
        Flowable.__init__(self)
        self.label = label
        self.bg = bg
        self.width = width
        self.height = height
        self.font_size = font_size
        self.text_color = text_color

    def wrap(self, availWidth, availHeight):
        return self.width, self.height

    def draw(self):
        c = self.canv
        r = self.height / 2.0
        c.saveState()
        c.setFillColor(self.bg)
        c.setStrokeColor(self.bg)
        c.roundRect(0, 0, self.width, self.height, r, stroke=1, fill=1)
        c.setFillColor(self.text_color)
        c.setFont("Helvetica-Bold", self.font_size)
        c.drawCentredString(self.width/2.0, self.height/2.0 - self.font_size*0.35, self.label[:32])
        c.restoreState()

class DrugVialIcon(Flowable):
    def __init__(self, label, solution_color=colors.HexColor("#4F46E5"), width=26*mm, height=42*mm, subtitle=""):
        Flowable.__init__(self)
        self.label = label
        self.solution_color = solution_color
        self.width = width
        self.height = height
        self.subtitle = subtitle

    def wrap(self, availWidth, availHeight):
        return self.width, self.height

    def draw(self):
        c = self.canv
        w = self.width
        h = self.height
        cx = w/2.0
        c.saveState()

        cap_h = h * 0.12
        neck_h = h * 0.08
        body_h = h * 0.56
        y0 = h * 0.18
        body_w = w * 0.64
        neck_w = body_w * 0.52

        c.setFillColor(colors.HexColor("#1F2937"))
        c.setStrokeColor(colors.HexColor("#1F2937"))
        c.roundRect(cx-body_w*0.36, y0+body_h+neck_h, body_w*0.72, cap_h, 1.2*mm, stroke=1, fill=1)

        glass = colors.HexColor("#E5F3FF")
        c.setFillColor(glass)
        c.setStrokeColor(colors.HexColor("#475569"))
        c.roundRect(cx-neck_w/2, y0+body_h, neck_w, neck_h, 1.2*mm, stroke=1, fill=1)
        c.roundRect(cx-body_w/2, y0, body_w, body_h, 2.4*mm, stroke=1, fill=1)

        liquid_h = body_h * 0.45
        c.setFillColor(self.solution_color)
        c.setStrokeColor(self.solution_color)
        c.roundRect(cx-body_w/2+1.0*mm, y0+1.0*mm, body_w-2.0*mm, liquid_h, 2.0*mm, stroke=0, fill=1)

        c.setFillColor(colors.white)
        c.rect(cx-body_w*0.26, y0+body_h*0.33, body_w*0.52, body_h*0.22, stroke=0, fill=1)
        c.setFillColor(colors.HexColor("#0F172A"))
        c.setFont("Helvetica-Bold", 7)
        c.drawCentredString(cx, y0+body_h*0.43, self.label[:10])

        c.setStrokeColor(colors.white)
        c.setLineWidth(1)
        c.line(cx-body_w*0.28, y0+body_h*0.10, cx-body_w*0.15, y0+body_h*0.18)
        c.line(cx-body_w*0.24, y0+body_h*0.18, cx-body_w*0.11, y0+body_h*0.26)

        c.setFillColor(colors.HexColor("#0F172A"))
        c.setFont("Helvetica-Bold", 8)
        c.drawCentredString(cx, 6, self.label[:16])
        if self.subtitle:
            c.setFont("Helvetica", 6.4)
            c.drawCentredString(cx, 0, self.subtitle[:18])

        c.restoreState()

def resistance_badge_color(category: str):
    cat = (category or "").strip().upper()
    if cat == "SENSITIVE STRAIN":
        return colors.HexColor("#16A34A")
    if cat in {"HR-TB", "RR-TB"}:
        return colors.HexColor("#D97706")
    if cat in {"MDR-TB", "PRE-XDR-TB", "XDR-TB"}:
        return colors.HexColor("#DC2626")
    return colors.HexColor("#475569")

def drug_color_map(drug: str):
    d = (drug or "").strip().lower()
    mapping = {
        "rifampicin": "#EF4444", "rifampin": "#EF4444", "isoniazid": "#3B82F6", "inh": "#3B82F6",
        "levofloxacin": "#8B5CF6", "moxifloxacin": "#7C3AED", "ofloxacin": "#A855F7", "gatifloxacin": "#9333EA",
        "bedaquiline": "#14B8A6", "linezolid": "#22C55E", "ethambutol": "#F59E0B", "pyrazinamide": "#EC4899",
        "amikacin": "#06B6D4", "kanamycin": "#0EA5E9", "capreomycin": "#6366F1",
    }
    return colors.HexColor(mapping.get(d, "#64748B"))

def make_summary_cards(sample_id, predicted, mixed, category):
    styles = getSampleStyleSheet()
    title_style = ParagraphStyle('TitleCard', parent=styles['Heading1'], fontName='Helvetica-Bold', fontSize=18, leading=22, textColor=colors.white, alignment=TA_LEFT)
    label_style = ParagraphStyle('LabelCard', parent=styles['Normal'], fontName='Helvetica', fontSize=8.5, leading=10, textColor=colors.white, alignment=TA_LEFT)
    value_style = ParagraphStyle('ValueCard', parent=styles['Normal'], fontName='Helvetica-Bold', fontSize=11.5, leading=13, textColor=colors.white, alignment=TA_LEFT)

    card_w = 55*mm

    left = Table([[Paragraph("Sample name", label_style)], [Paragraph(sample_id, title_style)]], colWidths=[card_w])
    left.setStyle(TableStyle([
        ('BACKGROUND', (0, 0), (-1, -1), colors.HexColor("#0F172A")),
        ('BOX', (0, 0), (-1, -1), 0.5, colors.HexColor("#0F172A")),
        ('LEFTPADDING', (0, 0), (-1, -1), 8), ('RIGHTPADDING', (0, 0), (-1, -1), 8),
        ('TOPPADDING', (0, 0), (-1, -1), 6), ('BOTTOMPADDING', (0, 0), (-1, -1), 7),
    ]))

    def card(title, value, bg):
        t = Table([[Paragraph(title, label_style)], [Paragraph(value, value_style)]], colWidths=[card_w])
        t.setStyle(TableStyle([
            ('BACKGROUND', (0, 0), (-1, -1), bg),
            ('BOX', (0, 0), (-1, -1), 0.5, bg),
            ('LEFTPADDING', (0, 0), (-1, -1), 8), ('RIGHTPADDING', (0, 0), (-1, -1), 8),
            ('TOPPADDING', (0, 0), (-1, -1), 6), ('BOTTOMPADDING', (0, 0), (-1, -1), 7),
        ]))
        return t

    return Table([[
        left,
        card("Predicted lineage", predicted, colors.HexColor("#1D4ED8")),
        card("Mixed lineage", mixed, colors.HexColor("#7C3AED")),
        card("Drug resistance", category, resistance_badge_color(category)),
    ]], colWidths=[card_w+1*mm]*4)

def build_lineage_table(lineage_references, results, matched_counts):
    styles = getSampleStyleSheet()
    hdr = ParagraphStyle('Hdr', parent=styles['Normal'], fontName='Helvetica-Bold', fontSize=9.2, textColor=colors.white, alignment=TA_CENTER)
    cell = ParagraphStyle('Cell', parent=styles['Normal'], fontName='Helvetica', fontSize=8.2, leading=9.5)
    lineages = sorted(lineage_references.keys())
    split = (len(lineages) + 1) // 2
    groups = [lineages[:split], lineages[split:]]

    def row_for_group(lin):
        if not lin:
            return ["", "", ""]
        return [Paragraph(lin, cell), Paragraph(f"{results.get(lin, 0.0):.2f}%", cell), Paragraph(matched_counts.get(lin, '0/0'), cell)]

    data = [[
        Paragraph("Lineage", hdr), Paragraph("Probability", hdr), Paragraph("Matched SNPs", hdr),
        Paragraph("Lineage", hdr), Paragraph("Probability", hdr), Paragraph("Matched SNPs", hdr),
    ]]
    max_rows = max(len(groups[0]), len(groups[1]))
    for i in range(max_rows):
        left = row_for_group(groups[0][i] if i < len(groups[0]) else None)
        right = row_for_group(groups[1][i] if i < len(groups[1]) else None)
        data.append(left + right)

    tbl = Table(data, colWidths=[46*mm, 24*mm, 24*mm, 46*mm, 24*mm, 24*mm], repeatRows=1)
    tbl.setStyle(TableStyle([
        ('BACKGROUND', (0, 0), (-1, 0), colors.HexColor("#0F766E")),
        ('TEXTCOLOR', (0, 0), (-1, 0), colors.white),
        ('GRID', (0, 0), (-1, -1), 0.25, colors.HexColor("#CBD5E1")),
        ('ROWBACKGROUNDS', (0, 1), (-1, -1), [colors.white, colors.HexColor("#F8FAFC")]),
        ('VALIGN', (0, 0), (-1, -1), 'MIDDLE'),
        ('LEFTPADDING', (0, 0), (-1, -1), 4), ('RIGHTPADDING', (0, 0), (-1, -1), 4),
        ('TOPPADDING', (0, 0), (-1, -1), 3), ('BOTTOMPADDING', (0, 0), (-1, -1), 3),
        ('LINEBEFORE', (3, 0), (3, -1), 0.9, colors.HexColor("#94A3B8")),
    ]))
    return tbl

def build_drug_hits_section(matched_res_csv: str, category: str):
    styles = getSampleStyleSheet()
    section_title = ParagraphStyle('SecTitle', parent=styles['Heading2'], fontName='Helvetica-Bold', fontSize=14, textColor=colors.HexColor("#0F172A"))
    body = [Paragraph("Drug-resistance overview", section_title), Spacer(1, 3*mm)]

    badge = Table([[PillBox(category, resistance_badge_color(category), width=34*mm, height=11*mm, font_size=8.4)]], colWidths=[36*mm])
    badge.setStyle(TableStyle([('LEFTPADDING', (0,0), (-1,-1), 0), ('RIGHTPADDING', (0,0), (-1,-1), 0)]))
    body.append(badge)
    body.append(Spacer(1, 2.5*mm))

    if not Path(matched_res_csv).exists():
        body.append(Paragraph("No resistance mutation file found.", styles['Normal']))
        return body

    df = pd.read_csv(matched_res_csv)
    if df.empty:
        note = Table([[Paragraph("Sensitive strain: no resistance-associated mutations detected.", ParagraphStyle('Sens', parent=styles['Normal'], fontName='Helvetica-Bold', fontSize=12, textColor=colors.white))]], colWidths=[120*mm])
        note.setStyle(TableStyle([
            ('BACKGROUND', (0, 0), (-1, -1), colors.HexColor("#16A34A")),
            ('LEFTPADDING', (0, 0), (-1, -1), 10), ('RIGHTPADDING', (0, 0), (-1, -1), 10),
            ('TOPPADDING', (0, 0), (-1, -1), 8), ('BOTTOMPADDING', (0, 0), (-1, -1), 8),
        ]))
        body.append(note)
        return body

    for c in ["drug", "Gene", "Mutation", "original_mutation", "confidence"]:
        if c not in df.columns:
            df[c] = ""
    df = df.copy()
    df["MutShown"] = df["Mutation"].fillna("")
    df.loc[df["MutShown"] == "", "MutShown"] = df["original_mutation"].fillna("")
    df = df.drop_duplicates(subset=["drug", "Gene", "MutShown", "confidence"])

    icon_cells = []
    detail_rows = [["Antibiotic", "Gene", "Mutation", "Prediction"]]
    for _, r in df.iterrows():
        drug = str(r.get("drug", "")).strip() or "Unknown"
        gene = str(r.get("Gene", "")).strip() or "-"
        mut = str(r.get("MutShown", "")).strip() or "-"
        conf = str(r.get("confidence", "")).strip() or "-"
        short = drug[:12]
        icon_cells.append(DrugVialIcon(short, solution_color=drug_color_map(drug), subtitle=gene))
        detail_rows.append([drug, gene, mut, conf])

    cols_per_row = 6
    icon_rows = [icon_cells[i:i+cols_per_row] for i in range(0, len(icon_cells), cols_per_row)]
    if icon_rows:
        last = icon_rows[-1]
        if len(last) < cols_per_row:
            last += [""] * (cols_per_row - len(last))
        flask_tbl = Table(icon_rows, colWidths=[28*mm] * cols_per_row)
        flask_tbl.setStyle(TableStyle([
            ('VALIGN', (0, 0), (-1, -1), 'TOP'), ('ALIGN', (0, 0), (-1, -1), 'CENTER'),
            ('LEFTPADDING', (0, 0), (-1, -1), 3), ('RIGHTPADDING', (0, 0), (-1, -1), 3),
            ('TOPPADDING', (0, 0), (-1, -1), 2), ('BOTTOMPADDING', (0, 0), (-1, -1), 3),
        ]))
        body.append(flask_tbl)
        body.append(Spacer(1, 3*mm))

    detail_tbl = Table(detail_rows, colWidths=[38*mm, 26*mm, 42*mm, 56*mm], repeatRows=1)
    detail_tbl.setStyle(TableStyle([
        ('BACKGROUND', (0, 0), (-1, 0), colors.HexColor("#7C3AED")),
        ('TEXTCOLOR', (0, 0), (-1, 0), colors.white),
        ('GRID', (0, 0), (-1, -1), 0.25, colors.HexColor("#CBD5E1")),
        ('ROWBACKGROUNDS', (0, 1), (-1, -1), [colors.white, colors.HexColor("#F8FAFC")]),
        ('VALIGN', (0, 0), (-1, -1), 'MIDDLE'),
        ('LEFTPADDING', (0, 0), (-1, -1), 5), ('RIGHTPADDING', (0, 0), (-1, -1), 5),
        ('TOPPADDING', (0, 0), (-1, -1), 4), ('BOTTOMPADDING', (0, 0), (-1, -1), 4),
        ('FONTSIZE', (0, 0), (-1, -1), 8.5),
        ('FONTNAME', (0, 0), (-1, 0), 'Helvetica-Bold'),
        ('FONTNAME', (0, 1), (-1, -1), 'Helvetica'),
    ]))
    body.append(detail_tbl)
    return body

def generate_infographic_pdf(sample_id: str, predicted: str, mixed: str, category: str,
                             lineage_references: dict, results: dict, matched_counts: dict,
                             matched_res_csv: str, output_pdf: str):
    doc = SimpleDocTemplate(
        output_pdf, pagesize=landscape(A3),
        leftMargin=10*mm, rightMargin=10*mm, topMargin=10*mm, bottomMargin=10*mm,
    )
    styles = getSampleStyleSheet()
    normal = ParagraphStyle('NormalReadable', parent=styles['Normal'], fontName='Helvetica', fontSize=10, leading=12)
    section_title = ParagraphStyle('SectionTitle', parent=styles['Heading2'], fontName='Helvetica-Bold', fontSize=14, textColor=colors.HexColor("#0F172A"))

    story = []
    story.append(make_summary_cards(sample_id, predicted, mixed, category))
    story.append(Spacer(1, 5*mm))
    story.append(Paragraph("Lineage profile", section_title))
    story.append(Spacer(1, 2*mm))
    story.append(KeepTogether(build_lineage_table(lineage_references, results, matched_counts)))
    story.append(Spacer(1, 5*mm))
    story.extend(build_drug_hits_section(matched_res_csv, category))
    story.append(Spacer(1, 4*mm))
    story.append(Paragraph("Readable infographic report generated in a compact single-page A3 layout for sample-level comparison.", normal))
    doc.build(story)
    logging.info(f"[{sample_id}] Infographic PDF saved to {output_pdf}")

def run_delly_call(sorted_bam: str, ref_genome: str, out_bcf: str, out_vcf: str, delly_cmd: str = "delly") -> bool:
    if not run_cmd(f"{delly_cmd} call -g {q(ref_genome)} -o {q(out_bcf)} {q(sorted_bam)}"):
        return False
    if not run_cmd(f"bcftools view {q(out_bcf)} -Ov -o {q(out_vcf)}"):
        return False
    return True

def merge_vcfs(vcf_inputs, out_vcf: str) -> bool:
    existing = [str(v) for v in vcf_inputs if v and Path(v).exists()]
    if not existing:
        logging.error("No input VCFs available to merge")
        return False
    if len(existing) == 1:
        if os.path.abspath(existing[0]) == os.path.abspath(out_vcf):
            return True
        return run_cmd(f"cp {q(existing[0])} {q(out_vcf)}")
    merged_list = " ".join(q(v) for v in existing)
    return run_cmd(f"bcftools concat -a -Ov -o {q(out_vcf)} {merged_list}")

def process_sample(sample_path, ref_genome, output_dir, lineage_references,
                   bed_file, threads=1, bam_override=None, vcf_override=None,
                   snpeff_cmd="snpEff", snpeff_db=DEFAULT_SNPEFF_DB,
                   mutations_csv=None, delly_cmd="delly", skip_delly=False,
                   make_infographic=True):
    start_time = time.time()
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    logs_dir = output_dir / "logs"
    logs_dir.mkdir(exist_ok=True)

    sample_id = Path(str(sample_path)).name
    sample_id = sample_id.split(",")[0]
    sample_id = Path(sample_id).stem
    sample_dir = output_dir / sample_id
    sample_dir.mkdir(exist_ok=True)

    vcf_lineage = None
    vcf_resist_gatk = None
    vcf_resist_delly = None
    vcf_resist_final = None
    sorted_bam = None
    final_category = "No data"

    if vcf_override and Path(vcf_override).exists():
        vcf_lineage = str(Path(vcf_override))
        vcf_resist_final = str(Path(vcf_override))
    else:
        if bam_override and Path(bam_override).exists():
            sorted_bam = str(Path(bam_override))
        else:
            kind, inputs = resolve_inputs(str(sample_path))
            if kind == "unknown":
                logging.error(f"[{sample_id}] Could not resolve inputs for '{sample_path}'.")
                return
            if kind == "vcf":
                vcf_lineage = inputs[0]
                vcf_resist_final = inputs[0]
            elif kind == "bam":
                sorted_bam = inputs[0]
            else:
                sam_file = sample_dir / f"{sample_id}.sam"
                bwa_log = logs_dir / f"{sample_id}.bwa.log"
                fastq_str = " ".join(q(x) for x in inputs)
                bwa_cmd = (
                    f"bwa mem -M -t {threads} -R '@RG\\tID:{sample_id}\\tSM:{sample_id}\\tPL:ILLUMINA' "
                    f"{q(ref_genome)} {fastq_str} > {q(sam_file)} 2> {q(bwa_log)}"
                )
                if not run_cmd(bwa_cmd) or not sam_file.exists():
                    return
                bam_file = sample_dir / f"{sample_id}.bam"
                sorted_bam_path = sample_dir / f"{sample_id}.sorted.bam"
                if not run_cmd(f"samtools view -S -b {q(sam_file)} -o {q(bam_file)}"):
                    return
                if not run_cmd(f"samtools sort -@ {threads} {q(bam_file)} -o {q(sorted_bam_path)}"):
                    return
                run_cmd(f"samtools index {q(sorted_bam_path)}")
                sorted_bam = str(sorted_bam_path)

        if sorted_bam and (not vcf_lineage or not vcf_resist_final):
            vcf_lineage = str(sample_dir / f"{sample_id}.lineage.vcf")
            tmp_bed = sample_dir / "targets.nc.bed"
            rewrite_bed_contig(bed_file, tmp_bed, target="NC_000962.3")
            if not run_cmd(f"gatk HaplotypeCaller -R {q(ref_genome)} -I {q(sorted_bam)} -O {q(vcf_lineage)} -L {q(tmp_bed)} --sample-ploidy 1"):
                return

            vcf_resist_gatk = str(sample_dir / f"{sample_id}.resist.gatk.vcf")
            if not run_cmd(f"gatk HaplotypeCaller -R {q(ref_genome)} -I {q(sorted_bam)} -O {q(vcf_resist_gatk)} --sample-ploidy 1"):
                return

            if not skip_delly:
                delly_bcf = str(sample_dir / f"{sample_id}.resist.delly.bcf")
                vcf_resist_delly = str(sample_dir / f"{sample_id}.resist.delly.vcf")
                if not run_delly_call(sorted_bam, ref_genome, delly_bcf, vcf_resist_delly, delly_cmd=delly_cmd):
                    vcf_resist_delly = None

            vcf_resist_final = str(sample_dir / f"{sample_id}.resist.combined.vcf")
            if not merge_vcfs([vcf_resist_gatk, vcf_resist_delly], vcf_resist_final):
                vcf_resist_final = vcf_resist_gatk

    sample_snps = set()
    with open(vcf_lineage, 'r') as vcf:
        for line in vcf:
            if line.startswith('#'):
                continue
            chrom, pos = line.strip().split('\t')[:2]
            chrom = contig_alias_to_nc(chrom)
            try:
                sample_snps.add((chrom, int(pos)))
            except ValueError:
                pass

    results, matched_counts = {}, {}
    for lineage, ref_snps in lineage_references.items():
        matched = len(ref_snps & sample_snps)
        total = len(ref_snps)
        results[lineage] = pct(matched, total)
        matched_counts[lineage] = f"{matched}/{total}"

    formatted = {lin: f"{results.get(lin, 0.0):.2f}%" for lin in lineage_references}
    sorted_lineages = sorted(results.items(), key=lambda x: x[1], reverse=True)

    if not sorted_lineages or sorted_lineages[0][1] == 0.0:
        predicted = "Unpredictable: insufficient SNP evidence"
    else:
        top = sorted_lineages[0][1]
        predicted_set = [lin for lin, pc in sorted_lineages if pc >= 90.0 or pc >= top - 25.0]
        predicted = ";".join(sorted(set(predicted_set)))

    high_conf = [lin for lin, pc in results.items() if pc >= 90.0]
    mixed = ";".join(sorted(high_conf)) if len(high_conf) > 1 else "-"

    elapsed_minutes = (time.time() - start_time) / 60
    lin_out = sample_dir / f"{sample_id}_lineage_result.txt"
    with open(lin_out, 'w') as f:
        f.write(f"Sample              : {sample_id}\n")
        f.write(f"Predicted lineage   : {predicted}\n")
        f.write(f"Mixed lineage       : {mixed}\n\n")
        f.write(f"{'Lineage':<20}{'Probability':<15}{'Matched SNPs'}\n")
        f.write("-" * 50 + "\n")
        for lineage in sorted(lineage_references.keys()):
            f.write(f"{lineage:<20}{formatted.get(lineage, '0.00%'):<15}{matched_counts.get(lineage, '0/0')}\n")
        f.write("\n" + "=" * 60 + "\n")
        f.write(f"Total Runtime: {elapsed_minutes:.2f} minutes\n")

    res_csv = str(sample_dir / "matched_mutations.csv")
    if mutations_csv and vcf_resist_final:
        ann_vcf = sample_dir / f"{sample_id}_ann.vcf"
        ann_csv = sample_dir / f"{sample_id}_annotated.csv"
        all_csv = sample_dir / "matched_all.csv"
        res_summary_csv = sample_dir / "mutations_summary.csv"
        final_txt = sample_dir / f"{sample_id}_final_prediction.txt"
        final_csv = sample_dir / f"{sample_id}_summary.csv"

        vcf_chr = sample_dir / f"{sample_id}.resist.chr.vcf"
        chrmap = sample_dir / "chr_rename.map"
        write_chrmap(str(chrmap), to_chr=True)
        if not run_cmd(f"bcftools annotate --rename-chrs {q(chrmap)} {q(vcf_resist_final)} -Ov -o {q(vcf_chr)}"):
            return
        if not run_cmd(f'{snpeff_cmd} {snpeff_db} {q(vcf_chr)} > {q(ann_vcf)}'):
            return

        df = parse_snpeff_vcf_to_df(str(ann_vcf))
        df.to_csv(ann_csv, index=False)
        match_catalogue(str(ann_csv), mutations_csv, str(all_csv), str(res_csv))
        write_drug_summary(str(res_csv), str(res_summary_csv))
        final_category = predict_final_category(str(res_csv), str(final_txt), str(final_csv))

    if make_infographic:
        pdf_out = sample_dir / f"{sample_id}_infographic_report.pdf"
        generate_infographic_pdf(
            sample_id=sample_id, predicted=predicted, mixed=mixed, category=final_category,
            lineage_references=lineage_references, results=results, matched_counts=matched_counts,
            matched_res_csv=str(res_csv), output_pdf=str(pdf_out),
        )

def main():
    ap = argparse.ArgumentParser(
        description="LineageXpress v8: A3 compact infographic PDF with improved drug icons and equal summary boxes",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    ap.add_argument("--sample_list", "--fastq", dest="sample_list", required=True, help="Text file: one sample spec per line ('R1,R2' | file/dir path | bare ID)")
    ap.add_argument("--ref_genome", default=DEFAULT_REF, help="Reference FASTA")
    ap.add_argument("--output_dir", default="results", help="Output directory")
    ap.add_argument("--snp_file", default=DEFAULT_SNP_FILE, help="Lineage SNP reference file")
    ap.add_argument("--bed_file", default=DEFAULT_BED_FILE, help="Target regions BED file for lineage VCF")
    ap.add_argument("--bam_file", default=None, help="Override: use this BAM for all samples")
    ap.add_argument("--vcf_file", default=None, help="Override: use this VCF for all samples")
    ap.add_argument("--threads", type=int, default=1, help="Threads for BWA / samtools")
    ap.add_argument("--mutations_csv", required=False, help="TB-Profiler WHO mutations.csv; enables mutation summary if provided")
    ap.add_argument("--snpeff_cmd", default="snpEff", help='Command for snpEff')
    ap.add_argument("--snpeff_db", default=DEFAULT_SNPEFF_DB, help="snpEff DB name")
    ap.add_argument("--delly_cmd", default="delly", help="Command for Delly executable")
    ap.add_argument("--skip_delly", action="store_true", help="Skip Delly structural variant calling")
    ap.add_argument("--no_infographic", action="store_true", help="Skip infographic PDF generation")
    args = ap.parse_args()

    for path, label in [(args.ref_genome, "ref_genome"), (args.snp_file, "snp_file"), (args.bed_file, "bed_file")]:
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
        process_sample(
            sample_path=sample, ref_genome=args.ref_genome, output_dir=args.output_dir,
            lineage_references=lineage_refs, bed_file=args.bed_file, threads=args.threads,
            bam_override=args.bam_file, vcf_override=args.vcf_file,
            snpeff_cmd=args.snpeff_cmd, snpeff_db=args.snpeff_db,
            mutations_csv=args.mutations_csv, delly_cmd=args.delly_cmd,
            skip_delly=args.skip_delly, make_infographic=(not args.no_infographic),
        )

if __name__ == "__main__":
    main()
