#!/usr/bin/env python3
"""
================================================================================
WP4-sc - MYELOID POPULATIONS IN REGENERATION vs SCARRING (H4)
GSE182141, single-cell, Acomys vs Mus, ear pinnae, d0/3/5/10/15
================================================================================

WHAT THIS REPLACES
------------------
H4 originally proposed INFERRING myeloid composition by deconvolving n=3 bulk
kidney samples per group, with the explicit caveat that deconvolution estimates
composition rather than observing it. GSE182141 observes it directly: 10X
single-cell, both species, five timepoints, from the group that defined the
regenerative phenotype (Adam, Potter & Seifert; PMID 38228141).

THE TISSUE CAVEAT - state this plainly in the manuscript
--------------------------------------------------------
This is EAR PINNAE, not kidney. Different organ, different injury model. The
4 mm ear punch is the canonical assay defining the deomyine regenerative
phenotype (Riddell et al. 2025 PNAS), so it is defensible - but H4 is then
answered in ear, and any extrapolation to kidney is inference, not result.

DATA FORMAT (verified 18 Aug 2026)
----------------------------------
CellRanger v2 layout, non-standard filenames:
    GSM5519169_Mus00_barcodes.tsv.gz   GSM5519174_Aco00_barcodes.tsv.gz
    GSM5519169_Mus00_genes.tsv.gz      GSM5519174_Aco00_genes.tsv.gz
    GSM5519169_Mus00_matrix.mtx.gz     GSM5519174_Aco00_matrix.mtx.gz

  * genes.tsv.gz is two columns: id<TAB>symbol
  * Mus 28,692 features; Acomys 24,509 - DIFFERENT references, so there is no
    shared feature space until symbols are joined.
  * Acomys IDs are mostly MOUSE SYMBOLS (Olfr1258, ...) from an Earlham
    Institute annotation; unnamed genes get ACOCA10068_EIv1_XXXXXXX. So a
    symbol intersection recovers most of the transcriptome.
  * ALL barcodes.tsv.gz are byte-identical (2,280,349 B) = the full 10X
    whitelist. These are RAW matrices. Cell calling and QC are on us.

Usage:
    python code\\04_scrna_myeloid.py --inspect      # report only, no analysis
    python code\\04_scrna_myeloid.py --load         # load + QC + save h5ad
    python code\\04_scrna_myeloid.py --score        # myeloid marker scoring

Requires: scanpy, anndata (pip install scanpy leidenalg igraph)
================================================================================
"""

from __future__ import annotations

import argparse
import gzip
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import config as cfg

log = cfg.log
GEO_DIR = cfg.RAW / "geo" / "GSE182141"
OUT = cfg.RESULTS / "wp4_scrna"

SAMPLES = [
    ("GSM5519169", "Mus00",   "Mus musculus",     0),
    ("GSM5519170", "Mus03",   "Mus musculus",     3),
    ("GSM5519171", "Mus05",   "Mus musculus",     5),
    ("GSM5519172", "Mus10",   "Mus musculus",    10),
    ("GSM5519173", "Mus15",   "Mus musculus",    15),
    ("GSM5519174", "Aco00",   "Acomys cahirinus", 0),
    ("GSM5519175", "Aco03",   "Acomys cahirinus", 3),
    ("GSM5519176", "Aco05",   "Acomys cahirinus", 5),
    ("GSM5519177", "Aco10",   "Acomys cahirinus",10),
    ("GSM5519178", "Aco15",   "Acomys cahirinus",15),
]

# ---- QC thresholds (pre-specified; raw matrices need cell calling) ----
MIN_COUNTS_PER_CELL = 500
MIN_GENES_PER_CELL = 200
MAX_PCT_MITO = 20.0
MIN_CELLS_PER_GENE = 3

# ---- marker panels ------------------------------------------------------------
# Deliberately literature-standard rather than tuned. Cross-species scoring is
# only as good as the symbols shared between the two references.
MARKERS = {
    # Panels deliberately over-specified: only genes present in BOTH references
    # are used, so each needs redundancy. The first inspect run showed Ptprc
    # (CD45) and Itgam (CD11b) absent from the shared space - the two markers
    # most datasets lean on - so myeloid identity rests on Csf1r/Adgre1/Fcgr1
    # plus the additions below.
    "myeloid_pan": ["Csf1r", "Fcgr1", "Adgre1", "Cd68", "Mertk", "Aif1",
                    "Cybb", "Ctss", "Tyrobp", "Fcer1g", "Lst1", "Lgals3"],
    "monocyte_classical": ["Ccr2", "Sell", "Plac8", "Ly6c1", "Vcan", "Fn1",
                           "S100a4", "Msrb1"],
    "monocyte_patrol": ["Cx3cr1", "Nr4a1", "Ace", "Spn", "Cd36", "Ear2",
                        "Pglyrp1"],
    "macro_proinflam": ["Il1b", "Tnf", "Nos2", "Cd86", "Socs3", "Cxcl2",
                        "Ccl3", "Ccl4", "Il6", "Nfkbia"],
    "macro_resolving": ["Mrc1", "Retnla", "Arg1", "Cd163", "Folr2", "Stab1",
                        "Il10", "Ccl24", "Cd209a", "Klf4"],
    "fibrogenic": ["Tgfb1", "Pdgfb", "Timp1", "Ctgf", "Serpine1", "Thbs1",
                   "Igf1"],
    "antifibrotic": ["Tgfb3", "Mmp9", "Mmp12", "Mmp13", "Mmp2", "Mmp14"],
    "fibroblast": ["Col1a1", "Col3a1", "Pdgfra", "Dcn", "Lum", "Postn",
                   "Fbn1"],
    "lymphoid": ["Cd3e", "Cd19", "Nkg7", "Ms4a1", "Cd8a", "Il7r", "Gzmb"],
    "endothelial": ["Pecam1", "Cdh5", "Kdr", "Egfl7"],
    "proliferating": ["Mki67", "Top2a", "Ccnb1", "Birc5"],
}

# A panel needs at least this many shared genes to be scored at all.
MIN_PANEL_GENES = 4


# ==============================================================================
# INSPECT
# ==============================================================================

def read_genes(path: Path) -> pd.DataFrame:
    rows = []
    with gzip.open(path, "rt") as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            rows.append((parts[0], parts[1] if len(parts) > 1 else parts[0]))
    return pd.DataFrame(rows, columns=["gene_id", "symbol"])


def inspect():
    log("=" * 74)
    log("GSE182141 INSPECTION - no analysis, just what is actually there")
    log("=" * 74)

    if not GEO_DIR.exists():
        log(f"Not found: {GEO_DIR}", "ERROR")
        log("  curl -L -o GSE182141_RAW.tar "
            '"https://www.ncbi.nlm.nih.gov/geo/download/'
            '?acc=GSE182141&format=file"', "ERROR")
        return False

    log(f"\n{'sample':10s} {'species':18s} {'day':>4s} {'features':>9s}  files")
    log("-" * 74)
    gene_tables = {}
    for gsm, name, species, day in SAMPLES:
        g = GEO_DIR / f"{gsm}_{name}_genes.tsv.gz"
        b = GEO_DIR / f"{gsm}_{name}_barcodes.tsv.gz"
        m = GEO_DIR / f"{gsm}_{name}_matrix.mtx.gz"
        have = "".join("+" if p.exists() else "-" for p in (b, g, m))
        n = "?"
        if g.exists():
            gt = read_genes(g)
            gene_tables[name] = gt
            n = f"{len(gt):,}"
        log(f"{name:10s} {species:18s} {day:>4d} {n:>9s}  {have}  (bar/gene/mtx)")

    mus = [v for k, v in gene_tables.items() if k.startswith("Mus")]
    aco = [v for k, v in gene_tables.items() if k.startswith("Aco")]
    if not (mus and aco):
        log("\nCould not read both species' gene tables.", "ERROR")
        return False

    m_sym = set(mus[0].symbol)
    a_sym = set(aco[0].symbol)
    shared = m_sym & a_sym

    log("\n" + "=" * 74)
    log("CROSS-SPECIES FEATURE SPACE")
    log("=" * 74)
    log(f"  Mus symbols          : {len(m_sym):,}")
    log(f"  Acomys symbols       : {len(a_sym):,}")
    log(f"  SHARED (join key)    : {len(shared):,}")
    log(f"  Acomys-only          : {len(a_sym - m_sym):,}")
    acoca = sum(1 for s in a_sym if s.startswith("ACOCA"))
    log(f"    of which ACOCA IDs : {acoca:,}  (unnamed in the EI annotation)")
    pct = 100 * len(shared) / max(len(a_sym) - acoca, 1)
    log(f"  -> {pct:.1f}% of NAMED Acomys genes have a Mus counterpart")

    log("\n" + "=" * 74)
    log("MARKER PANEL COVERAGE  (can each panel actually be scored?)")
    log("=" * 74)
    for panel, genes in MARKERS.items():
        both = [g for g in genes if g in shared]
        miss = [g for g in genes if g not in shared]
        flag = "" if len(both) >= MIN_PANEL_GENES else "   <-- TOO FEW TO SCORE"
        log(f"  {panel:20s} {len(both)}/{len(genes)} shared{flag}")
        if miss:
            log(f"      missing: {', '.join(miss)}")

    log("\n" + "=" * 74)
    log("CANONICAL MARKERS ABSENT FROM THE SHARED SPACE")
    log("=" * 74)
    for m in ("Ptprc", "Itgam", "Lyz2", "Itgax", "Chil3", "Spp1"):
        in_mus = m in m_sym
        in_aco = m in a_sym
        log(f"  {m:8s} Mus={'yes' if in_mus else 'NO ':3s}  "
            f"Acomys={'yes' if in_aco else 'NO ':3s}"
            f"{'   <-- absent from the Acomys annotation' if in_mus and not in_aco else ''}")
    log("  Genes present in Mus but not Acomys are annotation gaps, not")
    log("  biology. Do not interpret their absence as the cell type missing.")

    log("\n" + "=" * 74)
    log("FOCAL MODULE (from WP3) IN THIS DATASET")
    log("=" * 74)
    for g in list(cfg.FOCAL_GENES):
        # scRNA symbols are title-case; the focal set is upper-case
        cand = [g, g.capitalize(), g.title()]
        hit = next((c for c in cand if c in shared), None)
        log(f"  {g:8s} {'FOUND as ' + hit if hit else 'not in the shared space'}")

    log("\n" + "=" * 74)
    log("NOTE: barcodes files are the full 10X whitelist, so these are RAW")
    log("matrices. Cell calling and QC happen in --load, not inherited from")
    log("the authors. Expect most barcodes to be empty droplets.")
    log("=" * 74)
    return True


# ==============================================================================
# LOAD
# ==============================================================================

def load():
    try:
        import scanpy as sc
        import anndata as ad
        from scipy.io import mmread
        from scipy.sparse import csr_matrix
    except ImportError as e:
        log(f"missing package: {e}", "ERROR")
        log("  pip install scanpy leidenalg igraph", "ERROR")
        return

    OUT.mkdir(parents=True, exist_ok=True)
    adatas = []
    for gsm, name, species, day in SAMPLES:
        m = GEO_DIR / f"{gsm}_{name}_matrix.mtx.gz"
        g = GEO_DIR / f"{gsm}_{name}_genes.tsv.gz"
        b = GEO_DIR / f"{gsm}_{name}_barcodes.tsv.gz"
        if not all(p.exists() for p in (m, g, b)):
            log(f"  {name}: files missing, skipped", "WARN")
            continue

        log(f"  loading {name} ({species}, day {day})")
        with gzip.open(m, "rb") as fh:
            X = csr_matrix(mmread(fh).T)          # mtx is genes x cells
        genes = read_genes(g)
        with gzip.open(b, "rt") as fh:
            barcodes = [l.strip() for l in fh]

        a = ad.AnnData(X=X)
        a.var_names = pd.Index(genes.symbol.values).astype(str)
        a.var["gene_id"] = genes.gene_id.values
        a.obs_names = [f"{name}_{bc}" for bc in barcodes[:a.n_obs]]
        a.obs["sample"], a.obs["species"], a.obs["day"] = name, species, day
        a.var_names_make_unique()

        # ---- cell calling + QC (raw matrices) ----
        n0 = a.n_obs
        sc.pp.filter_cells(a, min_counts=MIN_COUNTS_PER_CELL)
        sc.pp.filter_cells(a, min_genes=MIN_GENES_PER_CELL)
        # Mitochondrial genes are named mt-* in the Mus reference but are
        # unnamed (ACOCA IDs) in the Acomys one. Filtering Mus on mito and
        # not Acomys applies DIFFERENT cell-quality thresholds per species,
        # which biases every downstream comparison. Record the fraction for
        # QC reporting, but do not filter on it.
        mito = [v for v in a.var_names if str(v).lower().startswith("mt-")]
        a.obs["pct_mito"] = (
            100 * np.asarray(a[:, mito].X.sum(1)).ravel()
            / np.maximum(np.asarray(a.X.sum(1)).ravel(), 1)
            if mito else np.nan)
        if not mito:
            log("    no mt- genes in this reference (expected for Acomys)")
        log(f"    {n0:,} barcodes -> {a.n_obs:,} cells")
        adatas.append(a)

    if not adatas:
        log("nothing loaded", "ERROR")
        return

    # ---- CRITICAL: restrict to genes shared by BOTH references -------------
    # An outer join pads species-specific genes with structural zeros. Since
    # the Mus reference names ~4,000 genes the Acomys one does not (including
    # Ptprc, Itgam, Lyz2, Spp1), per-cell totals then differ by species for
    # purely technical reasons - and normalize_total and score_genes both
    # depend on those totals. An inner join is the only defensible basis for
    # a cross-species comparison.
    shared = set(adatas[0].var_names)
    for a in adatas[1:]:
        shared &= set(a.var_names)
    shared = sorted(shared)
    log(f"\n  restricting to {len(shared):,} genes shared by all samples")
    log("  (inner join - species-specific genes are DROPPED, not zero-filled)")
    adatas = [a[:, shared].copy() for a in adatas]

    full = ad.concat(adatas, join="inner", label="batch",
                     keys=[a.obs["sample"].iloc[0] for a in adatas])
    sc.pp.filter_genes(full, min_cells=MIN_CELLS_PER_GENE)
    dest = OUT / "GSE182141_qc.h5ad"
    full.write_h5ad(dest)

    log("\n" + "=" * 74)
    log(f"  cells : {full.n_obs:,}")
    log(f"  genes : {full.n_vars:,}")
    log(f"  wrote : {dest}")
    log("\n  by sample:")
    for s, n in full.obs["sample"].value_counts().sort_index().items():
        log(f"    {s:8s} {n:>7,} cells")
    log("\n  NOTE: n = 1 animal per species per timepoint. Cells are NOT")
    log("  independent replicates. Any species difference here is descriptive;")
    log("  a formal test needs biological replicates that this dataset lacks.")
    log("=" * 74)
    log("  next: python code\\04_scrna_myeloid.py --score")


# ==============================================================================
# SCORE
# ==============================================================================

def score():
    try:
        import scanpy as sc
    except ImportError:
        log("pip install scanpy leidenalg igraph", "ERROR")
        return
    src = OUT / "GSE182141_qc.h5ad"
    if not src.exists():
        log(f"run --load first ({src} missing)", "ERROR")
        return

    a = sc.read_h5ad(src)
    log(f"loaded {a.n_obs:,} cells x {a.n_vars:,} genes")

    sc.pp.normalize_total(a, target_sum=1e4)
    sc.pp.log1p(a)

    log("\nscoring marker panels (shared symbols only):")
    for panel, genes in MARKERS.items():
        present = [g for g in genes if g in a.var_names]
        if len(present) < MIN_PANEL_GENES:
            log(f"  {panel:20s} SKIPPED - only {len(present)} genes present",
                "WARN")
            continue
        sc.tl.score_genes(a, present, score_name=panel)
        log(f"  {panel:20s} scored on {len(present)}/{len(genes)} genes")

    scored = [p for p in MARKERS if p in a.obs.columns]
    tab = (a.obs.groupby(["species", "day"], observed=True)[scored]
             .mean().round(3))
    OUT.mkdir(parents=True, exist_ok=True)
    tab.to_csv(OUT / "marker_scores_by_species_day.csv")
    a.write_h5ad(OUT / "GSE182141_scored.h5ad")

    log("\n" + "=" * 74)
    log("MEAN MARKER SCORE BY SPECIES AND DAY")
    log("=" * 74)
    log("\n" + tab.to_string())
    log(f"\nwrote {OUT / 'marker_scores_by_species_day.csv'}")
    log("\nThese are DESCRIPTIVE. Cross-species marker scores are sensitive to")
    log("reference and annotation differences between the two CellRanger runs.")
    log("Formal testing needs clustering, cell-type assignment, and a")
    log("proportion test with the animal as the unit - not the cell.")
    log("=" * 74)


def ratio():
    """Pre-registered secondary endpoint (PREREGISTRATION.md sec.3, H3):
    the CCR2:CX3CR1 trajectory.

    WHY A RATIO RATHER THAN A PANEL AVERAGE
    ---------------------------------------
    Ccr2 marks classical, fibrogenic monocytes; Cx3cr1 marks patrolling,
    non-fibrogenic ones. Their RATIO is directional and internally normalised:
    both genes come from the same cell, the same library, the same reference,
    so it is far less sensitive to the annotation and depth differences that
    make raw cross-species comparisons unreliable.

    Also computed:
      Tgfb1:Tgfb3 - the classical scarless-healing ratio, likewise
                    pre-registered and likewise internally controlled.

    STATISTICS. n = 1 animal per species per timepoint. Cells within a sample
    are pseudo-replicates, so a cell-level test would be spuriously
    significant. Per-sample summaries are reported with bootstrap CIs over
    CELLS, which quantify sampling precision WITHIN an animal and say nothing
    about between-animal variation. No p-value is produced, deliberately.
    """
    try:
        import scanpy as sc
    except ImportError:
        log("pip install scanpy", "ERROR")
        return
    src = OUT / "GSE182141_scored.h5ad"
    if not src.exists():
        src = OUT / "GSE182141_qc.h5ad"
    if not src.exists():
        log("run --load then --score first", "ERROR")
        return

    a = sc.read_h5ad(src)
    if "log1p" not in a.uns:
        sc.pp.normalize_total(a, target_sum=1e4)
        sc.pp.log1p(a)

    PAIRS = [("Ccr2", "Cx3cr1", "classical : patrolling monocyte"),
             ("Tgfb1", "Tgfb3", "fibrogenic : anti-scarring TGF-beta")]

    rng = np.random.default_rng(cfg.RANDOM_SEED)
    rows = []
    for num, den, label in PAIRS:
        if num not in a.var_names or den not in a.var_names:
            log(f"  {num}/{den}: not both present - skipped", "WARN")
            continue
        log(f"\n{'=' * 74}")
        log(f"{num} : {den}   ({label})")
        log("=" * 74)
        x = np.asarray(a[:, num].X.todense()).ravel() if hasattr(
            a[:, num].X, "todense") else np.asarray(a[:, num].X).ravel()
        y = np.asarray(a[:, den].X.todense()).ravel() if hasattr(
            a[:, den].X, "todense") else np.asarray(a[:, den].X).ravel()
        a.obs[f"_{num}"], a.obs[f"_{den}"] = x, y

        log(f"\n{'sample':8s} {'species':10s} {'day':>4s} {'cells':>7s} "
            f"{'%'+num:>8s} {'%'+den:>8s} {'log2 ratio':>11s}  95% CI (cells)")
        log("-" * 74)
        for s in sorted(a.obs["sample"].unique()):
            m = a.obs["sample"] == s
            xi, yi = x[m], y[m]
            n = int(m.sum())
            # mean expression per cell, pseudocount-stabilised
            r = np.log2((xi.mean() + 1e-3) / (yi.mean() + 1e-3))
            boot = [np.log2((xi[i].mean() + 1e-3) / (yi[i].mean() + 1e-3))
                    for i in (rng.integers(0, n, n) for _ in range(400))]
            lo, hi = np.percentile(boot, [2.5, 97.5])
            sp = "Acomys" if s.startswith("Aco") else "Mus"
            day = int(a.obs.loc[m, "day"].iloc[0])
            log(f"{s:8s} {sp:10s} {day:>4d} {n:>7,} "
                f"{100*(xi>0).mean():>7.1f}% {100*(yi>0).mean():>7.1f}% "
                f"{r:>11.3f}  [{lo:+.3f}, {hi:+.3f}]")
            rows.append({"pair": f"{num}:{den}", "sample": s, "species": sp,
                         "day": day, "n_cells": n,
                         f"pct_{num}_pos": 100*(xi>0).mean(),
                         f"pct_{den}_pos": 100*(yi>0).mean(),
                         "log2_ratio": r, "ci_lo": lo, "ci_hi": hi})

    if not rows:
        return
    df = pd.DataFrame(rows)
    OUT.mkdir(parents=True, exist_ok=True)
    df.to_csv(OUT / "focal_ratios_by_sample.csv", index=False)

    log("\n" + "=" * 74)
    log("INTERPRETATION GUARD")
    log("=" * 74)
    log("  The CIs above are bootstrapped over CELLS within one animal. They")
    log("  describe how precisely that animal was sampled - NOT whether the")
    log("  species differ. With n=1 animal per species per timepoint there is")
    log("  no between-animal variance to estimate, so no species comparison")
    log("  here can be tested. Non-overlapping CIs would NOT constitute")
    log("  evidence of a species difference.")
    log("")
    log(f"  wrote {OUT / 'focal_ratios_by_sample.csv'}")
    log("=" * 74)


def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--inspect", action="store_true")
    ap.add_argument("--load", action="store_true")
    ap.add_argument("--score", action="store_true")
    ap.add_argument("--ratio", action="store_true",
                    help="pre-registered CCR2:CX3CR1 trajectory (H3 secondary)")
    args = ap.parse_args()
    if args.inspect:
        inspect()
    elif args.load:
        load()
    elif args.score:
        score()
    elif args.ratio:
        ratio()
    else:
        ap.print_help()


if __name__ == "__main__":
    main()
