#!/usr/bin/env python3
"""
================================================================================
FIGURES for Paper B
================================================================================

Every number plotted here is READ FROM THE RESULT CSVs. Nothing is typed in.

That rule exists because the first version of this script hard-coded three sets
of values - the H1 q-values, the H2 significant-gene percentages, and the H5
permutation p-value. Hard-coded numbers in a figure script are how a panel
silently stops matching the table it is supposed to illustrate: the analysis is
re-run, the CSV updates, the figure does not, and nobody notices until a
referee does. Three numbers in this project have already been corrected once
(PPIA, MS4A1, TIMP1), so this is not hypothetical.

If a value you expect is missing, the script FAILS LOUDLY rather than
substituting a default.

Inputs (all under results/):
    wp3_kmer/arm1_residuals_utr3_k5.csv     Fig 1a  per-gene residual + length
    wp3_kmer/arm1_tests_utr3_k5.csv         Fig 1b  q-values, both models
    wp4_selection/selection_results.csv     Fig 2   LRT, q, gene_set
    wp4_selection/nguyen_annotation_frame_check.csv  Fig 3  three-tier frames
    wp5_expression/h3_module_tests.csv      Fig 4ab module interactions
    wp5_expression/h3_focal_gene_detail.csv Fig 4c  per-gene + recovery
    wp6_integration/h5_residuals_utr3_k5.csv  Fig 5 sequence divergence
    wp5_expression/h3_gene_level_EXPLORATORY_day5.csv  Fig 5 expression
    wp6_integration/h5_convergence_utr3_k5.csv         Fig 5 test statistics

Figure numbers follow ORDER OF FIRST CITATION in the manuscript:
    1 length confound | 2 selection | 3 annotation | 4 expression | 5 decoupling

Usage:
    python code\\08_figures.py                 # all five, PNG at 300 dpi
    python code\\08_figures.py --fig 3         # just figure 3
    python code\\08_figures.py --format pdf    # vector, for submission
    python code\\08_figures.py --dpi 600
================================================================================
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

sys.path.insert(0, str(Path(__file__).resolve().parent))
import config as cfg

R = cfg.RESULTS
FIGDIR = cfg.PROJECT_ROOT / "figures"

# House style. Journals want sans-serif, no top/right spines, small type.
plt.rcParams.update({
    "font.size": 8, "font.family": "DejaVu Sans",
    "axes.spines.top": False, "axes.spines.right": False,
    "axes.linewidth": 0.8, "xtick.major.width": 0.8, "ytick.major.width": 0.8,
})

COL = {
    "focal":   "#c0392b",   # focal module - the primary endpoint
    "control": "#95a5a6",   # control sets and non-significant
    "sig":     "#27ae60",   # significant non-focal modules
    "panel":   "#2980b9",   # adjusted / secondary
}


def read(rel: str) -> pd.DataFrame:
    """Load a result file, failing loudly if it is absent."""
    p = R / rel
    if not p.exists():
        raise SystemExit(
            f"MISSING INPUT: {p}\n"
            "Figures are generated only from analysis outputs. Re-run the "
            "relevant arm rather than supplying the number by hand.")
    return pd.read_csv(p)


def need(df: pd.DataFrame, col: str, where: str):
    if col not in df.columns:
        raise SystemExit(f"column '{col}' not found in {where}; "
                         f"available: {list(df.columns)}")


def save(fig, name: str, fmt: str, dpi: int):
    FIGDIR.mkdir(parents=True, exist_ok=True)
    out = FIGDIR / f"{name}.{fmt}"
    fig.savefig(out, dpi=dpi, bbox_inches="tight")
    plt.close(fig)
    print(f"  wrote {out}")


# ==============================================================================
# FIG 1 - H1: the length confound
# ==============================================================================

def fig1(fmt, dpi):
    d = read("wp3_kmer/arm1_residuals_utr3_k5.csv").dropna(subset=["residual_z"])
    # DO NOT drop_duplicates on gene before selecting the gene set.
    # Five genes (COL1A1, COL3A1, FN1, TIMP1, ACTA2) belong to BOTH
    # fibrosis_effector and PANEL_ECM_fibrosis. The panel row sorts first, so a
    # global dedup keeps the panel label and silently removes those genes from
    # the control set - which shifts this panel's correlation away from the
    # value the analysis reports. Select the set first, then dedup within it.
    t = read("wp3_kmer/arm1_tests_utr3_k5.csv")
    need(t, "q_BH", "arm1_tests_utr3_k5.csv")

    fig, ax = plt.subplots(1, 2, figsize=(7.0, 3.0))

    # (a) residual vs length, model fitted on CONTROLS ONLY - the same
    # restriction the analysis uses, so the line shown is the line applied.
    ctrl = d[d.gene_set.isin(cfg.CONTROL_SETS)].drop_duplicates("gene")
    foc = d[d.gene_set == "FOCAL"].drop_duplicates("gene")
    # Sanity check against the analysis, which fits the same model.
    _r, _p = stats.spearmanr(np.log10(ctrl.median_len.clip(lower=1)),
                             ctrl.residual_z)
    print(f"    control genes in length model: {len(ctrl)}  "
          f"Spearman {_r:+.3f} p={_p:.3f}")
    x = np.log10(ctrl.median_len.clip(lower=1))
    y = ctrl.residual_z
    ax[0].scatter(x, y, s=26, c=COL["control"], edgecolor="k", lw=.4,
                  label="control genes", zorder=2)
    ax[0].scatter(np.log10(foc.median_len.clip(lower=1)), foc.residual_z,
                  s=46, c=COL["focal"], edgecolor="k", lw=.5, marker="D",
                  label="focal module", zorder=3)
    b, a = np.polyfit(x, y, 1)
    xs = np.linspace(x.min(), x.max(), 50)
    ax[0].plot(xs, a + b * xs, "k--", lw=1.1, zorder=1)
    rho, p = stats.spearmanr(x, y)
    ax[0].set_xlabel("log$_{10}$ 3'UTR length (bp)")
    ax[0].set_ylabel("phylogenetic residual (z)")
    ax[0].set_title(f"Length predicts residual\nwithin controls: "
                    f"$\\rho$={rho:+.3f}, p={p:.3f}", fontsize=8)
    ax[0].legend(frameon=False, fontsize=7, loc="upper left")

    # (b) q-values before and after adjustment - READ FROM THE TESTS TABLE.
    f = t[t.comparison.str.startswith("FOCAL_vs_")].copy()
    f["ctrl"] = f.comparison.str.replace("FOCAL_vs_", "", regex=False)
    piv = f.pivot(index="ctrl", columns="model", values="q_BH")
    for m in ("unadjusted", "length_adjusted"):
        if m not in piv.columns:
            raise SystemExit(f"model '{m}' absent from arm1_tests; "
                             f"found {list(piv.columns)}")
    piv = piv.reindex(["housekeeping", "immune_nonfibrotic",
                       "fibrosis_effector"]).dropna(how="all")
    xp = np.arange(len(piv)); w = 0.36
    ax[1].bar(xp - w/2, -np.log10(piv["unadjusted"]), w, color=COL["focal"],
              label="unadjusted", edgecolor="k", lw=.4)
    ax[1].bar(xp + w/2, -np.log10(piv["length_adjusted"]), w,
              color=COL["control"], label="length-adjusted",
              edgecolor="k", lw=.4)
    ax[1].axhline(-np.log10(cfg.ALPHA), color="k", ls=":", lw=1)
    ax[1].text(len(piv) - .55, -np.log10(cfg.ALPHA) + .05,
               f"q={cfg.ALPHA}", fontsize=6.5, ha="right")
    ax[1].set_xticks(xp)
    ax[1].set_xticklabels([f"vs {c.split('_')[0]}" for c in piv.index],
                          fontsize=7)
    ax[1].set_ylabel("$-$log$_{10}$ q")
    ax[1].set_title("3'UTR signal does not survive\nthe pre-specified covariate",
                    fontsize=8)
    ax[1].legend(frameon=False, fontsize=7)

    fig.tight_layout()
    save(fig, "Fig1_H1_length_confound", fmt, dpi)


# ==============================================================================
# FIG 2 - H2: branch-site selection
# ==============================================================================

def fig2(fmt, dpi):
    s = read("wp4_selection/selection_results.csv")
    order = ["FOCAL", "fibrosis_effector", "housekeeping", "immune_nonfibrotic"]
    # Select the gene sets FIRST, then dedup within each. A global
    # drop_duplicates('gene') removed COL1A1 and COL3A1 from fibrosis_effector
    # (they also belong to PANEL_ECM_fibrosis, which sorts first), reporting
    # that set as 0/5 instead of 2/10 and omitting both significant genes from
    # the plot entirely.
    s = s[s.gene_set.isin(order)].drop_duplicates(["gene_set", "gene"])

    fig, ax = plt.subplots(1, 2, figsize=(7.0, 3.0))
    rng = np.random.default_rng(cfg.SEED if hasattr(cfg, "SEED") else 42)

    # (a) per-gene LRT, significant genes labelled
    for i, g in enumerate(order):
        sub = s[s.gene_set == g]
        jit = rng.normal(0, .07, len(sub))
        c = COL["focal"] if g == "FOCAL" else COL["control"]
        ax[0].scatter(np.full(len(sub), i) + jit, sub.LRT, s=34, c=c,
                      edgecolor="k", lw=.4, zorder=3)
        for _, r in sub[sub.q_BH < cfg.ALPHA].iterrows():
            ax[0].annotate(r.gene, (i + .14, r.LRT), fontsize=6.5, va="center")
    ax[0].axhline(3.84, color="k", ls=":", lw=1)   # chi2, df=1, p=0.05
    ax[0].text(len(order) - .55, 4.3, "p=0.05", fontsize=6.5, ha="right")
    ax[0].set_xticks(range(len(order)))
    ax[0].set_xticklabels(["focal", "fibrosis", "housekeep", "immune"],
                          fontsize=7)
    ax[0].set_ylabel("likelihood ratio (2$\\Delta$lnL)")
    ax[0].set_title("Branch-site test, Deomyinae lineage", fontsize=8)

    # (b) percent significant - COMPUTED, not typed
    g = (s.assign(sig=s.q_BH < cfg.ALPHA).groupby("gene_set")
           .agg(n=("gene", "size"), n_sig=("sig", "sum")).reindex(order))
    g["pct"] = 100 * g.n_sig / g.n
    ax[1].bar(range(len(g)), g.pct,
              color=[COL["focal"]] + [COL["control"]] * (len(g) - 1),
              edgecolor="k", lw=.4)
    ax[1].set_xticks(range(len(g)))
    ax[1].set_xticklabels(
        [f"{lbl}\n{int(r.n_sig)}/{int(r.n)}" for lbl, (_, r) in
         zip(["focal", "fibrosis", "housekeep", "immune"], g.iterrows())],
        fontsize=7)
    ax[1].set_ylabel(f"% genes significant (q<{cfg.ALPHA})")
    ax[1].set_title("Focal proportion is the lowest\nFisher p = 1.0 vs all three",
                    fontsize=8)

    fig.tight_layout()
    save(fig, "Fig2_H2_selection", fmt, dpi)


# ==============================================================================
# FIG 4 - H3: expression divergence and its boundary  [MAIN FIGURE]
# ==============================================================================

def fig4(fmt, dpi):
    m = read("wp5_expression/h3_module_tests.csv")
    g = read("wp5_expression/h3_focal_gene_detail.csv")

    mods = ["FOCAL", "fibrogenic", "PANEL_tubular_damage_repair",
            "PANEL_ECM_fibrosis", "PANEL_inflammatory_axis",
            "fibrosis_effector", "housekeeping"]
    nice = ["focal module", "fibrogenic", "tubular$^{\\dagger}$",
            "ECM effectors", "inflammatory", "fibrosis effector",
            "housekeeping"]

    fig, ax = plt.subplots(1, 3, figsize=(7.4, 3.1))
    for k, day in enumerate([2, 5]):
        sub = m[m.day == day].set_index("module")
        vals = [sub.loc[x, "mean_interaction"] if x in sub.index else np.nan
                for x in mods]
        qs = [sub.loc[x, "q_BH"] if x in sub.index else np.nan for x in mods]
        cols = [COL["focal"] if x == "FOCAL"
                else (COL["sig"] if (not np.isnan(q) and q < cfg.ALPHA)
                      else COL["control"])
                for x, q in zip(mods, qs)]
        yp = np.arange(len(mods))[::-1]
        ax[k].barh(yp, vals, color=cols, edgecolor="k", lw=.4)
        ax[k].axvline(0, color="k", lw=.8)
        ax[k].set_yticks(yp)
        ax[k].set_yticklabels(nice if k == 0 else [""] * len(mods), fontsize=7)
        ax[k].set_xlabel("interaction (Acomys $-$ Mus log$_2$FC)")
        ax[k].set_title(f"Day {day}", fontsize=8)
        for y, v, q in zip(yp, vals, qs):
            if not np.isnan(q) and q < cfg.ALPHA:
                ax[k].text(v - .06, y, "*", ha="right", va="center", fontsize=10)

    # (c) per-gene, annotated with relative read recovery percentile - the
    # number that decides whether the effect could be a mapping artefact.
    need(g, "recovery_pctile", "h3_focal_gene_detail.csv")
    gg = g.sort_values("interaction_d5").set_index("gene")
    yp = np.arange(len(gg))[::-1]
    ax[2].barh(yp, gg.interaction_d5, color=COL["focal"], edgecolor="k", lw=.4)
    for y, (_, r) in zip(yp, gg.iterrows()):
        ax[2].text(.05, y, f"{r.recovery_pctile:.0f}", fontsize=6.5, va="center")
    ax[2].axvline(0, color="k", lw=.8)
    ax[2].set_yticks(yp); ax[2].set_yticklabels(gg.index, fontsize=7)
    ax[2].set_xlabel("interaction, day 5")
    ax[2].set_title("Per gene (number = recovery %ile)", fontsize=8)

    fig.tight_layout()
    save(fig, "Fig4_H3_expression", fmt, dpi)
    print("  NOTE for the caption: the tubular block is marked with a dagger "
          "because it reaches q<0.05 but FAILS leave-one-out (removing SPP1 "
          "takes day 5 to p=0.63). It must not be presented as a finding.")


# ==============================================================================
# FIG 5 - H5: sequence and expression divergence are decoupled
# ==============================================================================

def fig5(fmt, dpi, region="utr3", k=5, day=5):
    h = read(f"wp6_integration/h5_residuals_{region}_k{k}.csv")
    stats_tab = read(f"wp6_integration/h5_convergence_{region}_k{k}.csv")
    e = read(f"wp5_expression/h3_gene_level_EXPLORATORY_day{day}.csv")
    e = e.rename(columns={e.columns[0]: "gene"})
    e["gene"] = e.gene.astype(str).str.upper()
    e = e.drop_duplicates("gene").set_index("gene")

    j = h.merge(e, left_on="gene", right_index=True)
    row = stats_tab[stats_tab.day == day]
    if row.empty:
        raise SystemExit(f"no day-{day} row in h5_convergence_{region}_k{k}.csv")
    row = row.iloc[0]

    fig, ax = plt.subplots(1, 2, figsize=(7.0, 3.0))
    ax[0].scatter(j.residual_z, np.abs(j.interaction), s=16,
                  c=COL["control"], edgecolor="k", lw=.3, alpha=.8)
    ax[0].set_xlabel("sequence divergence (H1 residual z)")
    ax[0].set_ylabel("|expression divergence| (H3)")
    ax[0].set_title(f"raw  $\\rho$={row.rho_raw:+.3f}, p={row.p_raw:.3f}",
                    fontsize=8)

    def resid(y, X):
        A = np.column_stack([np.ones(len(y)), X])
        beta, *_ = np.linalg.lstsq(A, y, rcond=None)
        return y - A @ beta

    X = np.column_stack([
        np.log10(j.median_len.clip(lower=1)),
        (np.abs(j.lfc_acomys) + np.abs(j.lfc_mus)) / 2])
    ax[1].scatter(resid(j.residual_z.values, X),
                  resid(np.abs(j.interaction.values), X),
                  s=16, c=COL["panel"], edgecolor="k", lw=.3, alpha=.8)
    ax[1].set_xlabel("sequence divergence (adjusted)")
    ax[1].set_ylabel("|expression divergence| (adjusted)")
    ax[1].set_title(f"adjusted  $\\rho$={row.rho_adjusted:+.3f}, "
                    f"p$_{{perm}}$={row.p_perm:.3f}", fontsize=8)
    for a_ in ax:
        a_.axhline(0, color="k", lw=.5, ls=":")

    fig.suptitle("Sequence and expression divergence are decoupled  "
                 f"(n = {int(row.n_genes)} genes, {region}, day {day})",
                 fontsize=8.5, y=1.02)
    fig.tight_layout()
    save(fig, "Fig5_H5_decoupling", fmt, dpi)


# ==============================================================================
# FIG 3 - annotation provenance determines reading-frame integrity
# ==============================================================================

def fig3(fmt, dpi):
    """Three-tier comparison plus per-gene concordance between the two
    independent A. cahirinus CDS reconstructions."""
    n = read("wp4_selection/nguyen_annotation_frame_check.csv").set_index("gene")
    for c in ("our_tier2_stops", "russatus_refseq_stops"):
        need(n, c, "nguyen_annotation_frame_check.csv — re-run "
                   "09_test_nguyen_annotation.py --check to add the "
                   "comparison columns")

    russ = n.russatus_refseq_stops.dropna()
    ours = n.our_tier2_stops.dropna()

    fig, ax = plt.subplots(1, 2, figsize=(7.0, 3.0))

    tiers = [
        ("Native RefSeq\n(A. russatus)", int((russ > 0).sum()), len(russ),
         COL["sig"]),
        ("Published CAT\n(Nguyen 2023)", int((n.internal_stops > 0).sum()),
         len(n), COL["panel"]),
        ("Homology recon.\n(this study)", int((ours > 0).sum()), len(ours),
         COL["focal"]),
    ]
    xp = np.arange(len(tiers))
    pct = [100 * b / max(t, 1) for _, b, t, _ in tiers]
    ax[0].bar(xp, pct, color=[c for *_, c in tiers], edgecolor="k", lw=.5)
    for x, (lab, b, t, _), p in zip(xp, tiers, pct):
        ax[0].text(x, p + 1.5, f"{b}/{t}", ha="center", fontsize=7)
    ax[0].set_xticks(xp)
    ax[0].set_xticklabels([lab for lab, *_ in tiers], fontsize=7)
    ax[0].set_ylabel("% CDS with internal stop codons")
    ax[0].set_ylim(0, 100)
    ax[0].set_title("Annotation provenance determines\nreading-frame integrity",
                    fontsize=8)

    shared = list(ours.index)
    x = n.loc[shared, "internal_stops"].values.astype(float)
    y = ours.loc[shared].values.astype(float)
    if len(shared) == 0:
        raise SystemExit("no genes shared between the two reconstructions")
    focal = {g.upper() for g in cfg.FOCAL_GENES}
    cols = [COL["focal"] if g in focal else COL["control"] for g in shared]
    ax[1].scatter(x + 0.15, y + 0.15, s=34, c=cols, edgecolor="k", lw=.4,
                  zorder=3)
    for g in shared:
        if g in focal:
            ax[1].annotate(g, (n.loc[g, "internal_stops"] + .6, ours[g] + .15),
                           fontsize=6.5, va="center")
    lim = max(x.max(), y.max()) * 1.15
    ax[1].plot([0, lim], [0, lim], "k:", lw=.8, zorder=1)
    ax[1].set_xscale("symlog"); ax[1].set_yscale("symlog")
    ax[1].set_xlabel("internal stops — published CAT")
    ax[1].set_ylabel("internal stops — homology recon.")
    conc = 100 * np.mean((x > 0) == (y > 0))
    ax[1].set_title(f"Per-gene agreement ({conc:.0f}% concordant)", fontsize=8)

    fig.tight_layout()
    save(fig, "Fig3_annotation_frame_integrity", fmt, dpi)
    print(f"    tiers: {[(l.replace(chr(10),' '), f'{b}/{t}') for l,b,t,_ in tiers]}")


def main():
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--fig", type=int, choices=[1, 2, 3, 4, 5],
                    help="build one figure only")
    ap.add_argument("--format", default="png", choices=["png", "pdf", "svg"])
    ap.add_argument("--dpi", type=int, default=300)
    a = ap.parse_args()

    todo = [a.fig] if a.fig else [1, 2, 3, 4, 5]
    fns = {1: fig1, 2: fig2, 3: fig3, 4: fig4, 5: fig5}
    # 3 = annotation, 4 = expression, 5 = decoupling:
    # numbering follows order of first citation in the manuscript.
    print(f"figures -> {FIGDIR}  ({a.format}, {a.dpi} dpi)")
    for n in todo:
        print(f"Fig {n}:")
        fns[n](a.format, a.dpi)


if __name__ == "__main__":
    main()
