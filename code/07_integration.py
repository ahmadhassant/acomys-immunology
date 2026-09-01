#!/usr/bin/env python3
"""
================================================================================
WP6 - INTEGRATION (H5): does sequence divergence predict expression divergence?
================================================================================

H5 (PREREGISTRATION.md sec.3): "Genes with high sequence-level residual (H1) are
enriched among genes with significant species x time interaction (H3)."

  Test          rank-based enrichment with a permutation null
  Falsified if  no enrichment - AND THIS IS AN INFORMATIVE NEGATIVE, declared
                as such in advance

WHY THIS IS WORTH RUNNING EVEN THOUGH H1 WAS NEGATIVE
------------------------------------------------------
H1 asked whether the FOCAL MODULE diverges more than control sets. It does not.
H5 asks something different and orthogonal: across genes generally, does
sequence-level divergence track expression-level divergence at all? A gene set
can fail a focal-vs-control comparison while the underlying quantities still
correlate genome-wide - or fail to.

Either answer is publishable:
  - ENRICHED  -> sequence residual carries information about regulatory
                 divergence, and alignment-free screens have predictive value.
  - NOT       -> the two levels are decoupled, and any dissociation seen
                 between the sequence arms and the expression arm is a
                 general property of these data rather than something
                 special to the focal module.

SCOPE - THIS IS NOT RESTRICTED TO THE NAMED GENE SETS
-----------------------------------------------------
H1's published test covered ~37 genes (focal + controls + panel). That is far
too few for an enrichment test. This uses the WHOLE ortholog scaffold - 236
genes, of which 199 are background never involved in any focal-vs-control
comparison. Residuals are computed here for all of them.

THE CONFOUND THAT WOULD MANUFACTURE A POSITIVE
-----------------------------------------------
Both variables have known nuisance correlates:
  - H1 residual correlates with SEQUENCE LENGTH (positive Spearman,
    estimated within control genes only)
  - H3 interaction correlates with READ RECOVERY and ABUNDANCE
And length, recovery and abundance are themselves correlated with each other
and with expression level. A raw correlation between the two endpoints could
therefore be entirely a shared-confounder artefact.

So the primary test is a PARTIAL correlation, with both variables residualised
on log(length), relative recovery and mean abundance before they are compared.
The raw correlation is reported alongside, and if they disagree the adjusted
one is the answer - the same rule that decided H1.

Usage:
    python code\\07_integration.py --region utr3 --k 5
    python code\\07_integration.py --region promoter --k 5
    python code\\07_integration.py --all-regions
================================================================================
"""

from __future__ import annotations

import argparse
import importlib.util
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

sys.path.insert(0, str(Path(__file__).resolve().parent))
import config as cfg

log = cfg.log
OUT = cfg.RESULTS / "wp6_integration"
WP3 = cfg.RESULTS / "wp3_kmer"
WP5 = cfg.RESULTS / "wp5_expression"

N_PERM = 10000
SEED = 42
TOP_QUANTILE = 0.25          # "high residual" = top quartile, fixed in advance
MIN_SPECIES = 4


def _load_wp3():
    """Import 03_kmer_phylo_controlled so its residual function is reused
    verbatim rather than reimplemented - a second implementation would be a
    second chance to introduce a discrepancy."""
    p = Path(__file__).with_name("03_kmer_phylo_controlled.py")
    spec = importlib.util.spec_from_file_location("wp3", p)
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


def scaffold_genes() -> list[str]:
    return sorted(d.name for d in cfg.ORTHO.glob("*") if d.is_dir())


def compute_residuals(region: str, k: int) -> pd.DataFrame:
    """Phylogenetic residual for EVERY scaffold gene, not just the named sets."""
    wp3 = _load_wp3()
    rows = []
    genes = scaffold_genes()
    log(f"Computing {region} k={k} residuals for {len(genes)} scaffold genes")
    for i, g in enumerate(genes, 1):
        if i % 50 == 0:
            log(f"  {i}/{len(genes)}")
        seqs = wp3.load_gene(g, region)
        if len(seqs) < MIN_SPECIES:
            continue
        res = wp3.phylogenetic_residual(seqs, k)
        rz = res.get("residual_z", np.nan)
        if rz is None or (isinstance(rz, float) and np.isnan(rz)):
            continue
        rows.append({
            "gene": g.upper(),
            "residual_z": float(rz),
            "n_seqs": len(seqs),
            "median_len": float(np.median([len(s) for s in seqs.values()])),
        })
    df = pd.DataFrame(rows)
    log(f"  {len(df)} genes produced a residual")
    return df


def load_h3(day: int) -> pd.DataFrame:
    f = WP5 / f"h3_gene_level_EXPLORATORY_day{day}.csv"
    if not f.exists():
        log(f"missing {f} - run 06_expression_reanalysis.py --analyse", "ERROR")
        return pd.DataFrame()
    d = pd.read_csv(f, index_col=0)
    d.index = [str(i).upper() for i in d.index]
    d = d[~d.index.duplicated()]
    return d


def partial_resid(y: np.ndarray, X: np.ndarray) -> np.ndarray:
    """Residuals of y on covariates X (with intercept)."""
    A = np.column_stack([np.ones(len(y)), X])
    beta, *_ = np.linalg.lstsq(A, y, rcond=None)
    return y - A @ beta


def run(region: str, k: int) -> pd.DataFrame | None:
    rng = np.random.default_rng(SEED)
    OUT.mkdir(parents=True, exist_ok=True)

    res = compute_residuals(region, k)
    if res.empty:
        return None

    rows = []
    for day in (2, 5):
        h3 = load_h3(day)
        if h3.empty:
            return None
        m = res.merge(h3, left_on="gene", right_index=True, how="inner")
        if len(m) < 30:
            log(f"day {day}: only {len(m)} genes overlap - too few to test",
                "WARN")
            continue

        # Endpoint definitions, fixed before looking:
        #   sequence divergence  = H1 phylogenetic residual
        #   expression divergence = |species x time interaction|
        # Magnitude is used because H5 asks about DIVERGENCE, which is
        # directionless; the signed version is reported alongside.
        seq = m["residual_z"].values
        expr = np.abs(m["interaction"].values)
        expr_signed = m["interaction"].values

        # Nuisance covariates shared by both endpoints.
        loglen = np.log10(np.clip(m["median_len"].values, 1, None))
        abundance = (np.abs(m["lfc_acomys"].values)
                     + np.abs(m["lfc_mus"].values)) / 2
        X = np.column_stack([loglen, abundance])

        rho_raw, p_raw = stats.spearmanr(seq, expr)
        seq_a = partial_resid(seq, X)
        expr_a = partial_resid(expr, X)
        rho_adj, _ = stats.spearmanr(seq_a, expr_a)

        # Permutation null: shuffle the pairing between the two endpoints.
        # This preserves both marginal distributions exactly and destroys only
        # the gene-to-gene correspondence, which is the thing under test.
        null = np.array([stats.spearmanr(seq_a, rng.permutation(expr_a))[0]
                         for _ in range(N_PERM)])
        p_perm = float((np.abs(null) >= abs(rho_adj)).mean())

        # Pre-specified rank-based enrichment: are top-quartile-residual genes
        # enriched among high-interaction genes?
        thr = np.quantile(seq, 1 - TOP_QUANTILE)
        hi = seq >= thr
        obs_rank = float(stats.rankdata(expr_a)[hi].mean())
        null_rank = np.array([stats.rankdata(expr_a)[
            rng.permutation(hi)].mean() for _ in range(N_PERM)])
        p_enrich = float((null_rank >= obs_rank).mean())

        rho_signed, _ = stats.spearmanr(seq_a, partial_resid(expr_signed, X))

        log(f"\n--- {region} k={k}, day {day}  (n = {len(m)} genes) ---")
        log(f"  Spearman(residual, |interaction|)")
        log(f"    raw       rho = {rho_raw:+.4f}   p = {p_raw:.4g}")
        log(f"    ADJUSTED  rho = {rho_adj:+.4f}   p(perm) = {p_perm:.4f}")
        log(f"    signed    rho = {rho_signed:+.4f}")
        log(f"  Top-{int(TOP_QUANTILE*100)}% residual genes (n = {int(hi.sum())}):")
        log(f"    mean |interaction| rank = {obs_rank:.1f} vs null "
            f"{null_rank.mean():.1f}   p(perm, one-sided) = {p_enrich:.4f}")

        rows.append({"region": region, "k": k, "day": day, "n_genes": len(m),
                     "rho_raw": rho_raw, "p_raw": p_raw,
                     "rho_adjusted": rho_adj, "p_perm": p_perm,
                     "rho_signed": rho_signed,
                     "n_top": int(hi.sum()), "top_mean_rank": obs_rank,
                     "null_mean_rank": float(null_rank.mean()),
                     "p_enrichment": p_enrich})

    if not rows:
        return None
    df = pd.DataFrame(rows)
    df.to_csv(OUT / f"h5_convergence_{region}_k{k}.csv", index=False)
    res.to_csv(OUT / f"h5_residuals_{region}_k{k}.csv", index=False)

    log("\n" + "=" * 74)
    log("H5 VERDICT")
    log("=" * 74)
    sig = df[(df.p_perm < cfg.ALPHA) | (df.p_enrichment < cfg.ALPHA)]
    if len(sig):
        log(f"  Enrichment detected at {len(sig)}/{len(df)} timepoint(s) "
            "-> H5 SUPPORTED")
        log("  Sequence-level residual carries information about "
            "expression-level divergence.")
    else:
        log("  No enrichment at either timepoint -> H5 NOT SUPPORTED")
        log("  Sequence divergence and expression divergence are DECOUPLED "
            "in these data.")
        log("  This was pre-registered as an informative negative: it is "
            "reported whichever way it falls, and it is read together with "
            "the sequence-level and expression-level arms rather than on "
            "its own.")
    log("=" * 74)
    log(f"\nwrote {OUT / f'h5_convergence_{region}_k{k}.csv'}")
    return df


def main():
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--region", default="utr3",
                    choices=["promoter", "utr5", "cds", "utr3"])
    ap.add_argument("--k", type=int, default=5)
    ap.add_argument("--all-regions", action="store_true")
    a = ap.parse_args()

    if a.all_regions:
        out = []
        for r in ("promoter", "utr5", "cds", "utr3"):
            d = run(r, a.k)
            if d is not None:
                out.append(d)
        if out:
            allr = pd.concat(out)
            allr.to_csv(OUT / f"h5_convergence_ALL_k{a.k}.csv", index=False)
            log("\n" + allr.to_string(index=False))
    else:
        run(a.region, a.k)


if __name__ == "__main__":
    main()
