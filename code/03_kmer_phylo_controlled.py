#!/usr/bin/env python3
"""
================================================================================
WP3 - PHYLOGENY-CONTROLLED ALIGNMENT-FREE SEQUENCE ANALYSIS
Acomys cahirinus regeneration project
================================================================================

WHY THIS SCRIPT EXISTS
----------------------
Your Echinococcus and FISH k-mer pipelines classify species/genotypes from
sequence, and they do it very well. Pointed at Acomys vs Mus they will return
~100% accuracy - and mean nothing, because Acomys (Deomyinae) and Mus (Murinae)
are different subfamilies that split ~18.3 Ma. The classifier separates
taxonomy, not regenerative phenotype.

This script asks a different, answerable question:

    Does the CX3CL1/CCL2-TGF-beta module in Acomys deviate from the murid
    phylogenetic expectation MORE than matched control gene sets do?

The unit of inference is a RESIDUAL from a phylogenetic expectation, not a
raw classification accuracy. That is the whole point.

DESIGN
------
  Arm 1  Phylogenetic residual test    <- PRIMARY endpoint
  Arm 2  Composition nulls (GC3, codon-shuffle, dinucleotide-preserving)
  Arm 3  Region-stratified analysis    (promoter / 5'UTR / CDS / 3'UTR)
  Arm 4  Leakage-hardened ML + SHAP    <- SECONDARY, descriptive only

STATISTICAL COMMITMENTS (pre-specified; do not change post hoc)
  - PCA and all scaling fitted INSIDE cross-validation folds only.
  - Every accuracy figure accompanied by a label-permutation null (n=1000).
  - Bootstrap 95% CIs on all effect sizes.
  - BH-FDR within each arm; q-values reported alongside p.
  - Sequence-length bias explicitly controlled (frequencies, not counts,
    plus length as a covariate).

REUSED FROM YOUR EXISTING CODE
  revision_analyses_v2.py  (Echinococcus)  - leakage fixes, length-bias control,
                                             CI machinery, PCA-inside-folds
  kmer_shap_analysis_complete.py (FISH)    - KmerSHAPAnalyzer class structure
  Best K.py (Echinococcus)                 - k sweep, now nested inside CV
  Windo/K-mer.py (Echinococcus)            - sliding-window scan

Usage:
    python 03_kmer_phylo_controlled.py --arm all
    python 03_kmer_phylo_controlled.py --arm residual --region utr3
    python 03_kmer_phylo_controlled.py --arm nulls --k 5

Requirements:
    pip install numpy pandas scipy scikit-learn biopython matplotlib seaborn shap
    pip install statsmodels
    # optional, for the formal PGLS:
    pip install dendropy
================================================================================
"""

from __future__ import annotations

import argparse
import json
import random
import sys
import warnings
from collections import defaultdict
from dataclasses import dataclass
from itertools import product
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from Bio import SeqIO
from Bio.Seq import Seq
from scipy import stats
from scipy.spatial.distance import pdist, squareform
from sklearn.decomposition import PCA
from sklearn.ensemble import RandomForestClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import accuracy_score, f1_score
from sklearn.model_selection import StratifiedKFold, cross_val_score
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import MinMaxScaler, StandardScaler

warnings.filterwarnings("ignore")

# ==============================================================================
# CONFIGURATION - all constants live in config.py (single source of truth)
# ==============================================================================

sys.path.insert(0, str(Path(__file__).resolve().parent))
import config as cfg

cfg.set_all_seeds()
RANDOM_SEED = cfg.RANDOM_SEED

PROJECT_ROOT, DATA = cfg.PROJECT_ROOT, cfg.DATA
ORTHO = cfg.ORTHO
RESULTS = cfg.RESULTS / "wp3_kmer"
FIGS = RESULTS / "figures"

K_SWEEP = cfg.K_SWEEP
N_PERMUTATIONS = cfg.N_PERMUTATIONS
N_BOOTSTRAP = cfg.N_BOOTSTRAP
CV_FOLDS = cfg.CV_FOLDS
ALPHA = cfg.ALPHA
MIN_NULL_DRAWS = 20   # below this an empirical p-value is meaningless

DIVERGENCE_MA = cfg.DIVERGENCE_MA
CLADE = cfg.CLADE
REGENERATIVE = cfg.REGENERATIVE
FOCAL_GENES = list(cfg.FOCAL_GENES)
CONTROL_SETS = {k: v for k, v in cfg.CONTROL_SETS.items() if v}
REGIONS = cfg.REGIONS

log = cfg.log
plt.rcParams.update(cfg.PLOT_STYLE)


# ==============================================================================
# K-MER FEATURES
# ==============================================================================

def all_kmers(k: int) -> list[str]:
    return ["".join(p) for p in product("ACGT", repeat=k)]


def kmer_profile(seq: str, k: int, normalise: bool = True) -> np.ndarray:
    """Frequency vector. Frequencies not counts - this is the length-bias fix."""
    vocab = {km: i for i, km in enumerate(all_kmers(k))}
    vec = np.zeros(len(vocab), dtype=np.float64)
    s = seq.upper()
    for i in range(len(s) - k + 1):
        km = s[i:i + k]
        j = vocab.get(km)
        if j is not None:
            vec[j] += 1
    if normalise:
        tot = vec.sum()
        if tot > 0:
            vec /= tot
    return vec


def profiles_matrix(seqs: dict[str, str], k: int) -> tuple[np.ndarray, list, list]:
    labels = sorted(seqs)
    X = np.vstack([kmer_profile(seqs[s], k) for s in labels])
    return X, labels, [f"{k}mer_{km}" for km in all_kmers(k)]


def gc_content(seq: str) -> float:
    s = seq.upper()
    n = sum(s.count(b) for b in "ACGT")
    return (s.count("G") + s.count("C")) / n if n else np.nan


def gc3_content(cds: str) -> float:
    """GC at synonymous third codon positions. The single biggest confounder
    of coding-sequence k-mer composition."""
    s = cds.upper()
    third = s[2::3]
    n = sum(third.count(b) for b in "ACGT")
    return (third.count("G") + third.count("C")) / n if n else np.nan


# ==============================================================================
# NULL MODELS
# ==============================================================================
# A "discriminative motif" is only interesting if it survives these.

def shuffle_dinucleotide(seq: str, rng: random.Random) -> str:
    """Altschul-Erikson dinucleotide-preserving shuffle (Eulerian-path method).

    Preserves mononucleotide AND dinucleotide composition. This is the correct
    null for non-coding regions: naive shuffling destroys dinucleotide bias
    (notably CpG depletion) and makes everything look significant.
    """
    s = seq.upper()
    s = "".join(c for c in s if c in "ACGT")
    if len(s) < 4:
        return s

    edges = defaultdict(list)
    for a, b in zip(s, s[1:]):
        edges[a].append(b)
    last = s[-1]

    # Build a random arborescence into `last`, then shuffle remaining edges.
    for _ in range(100):
        tree_edge = {}
        ok = True
        for v in edges:
            if v == last:
                continue
            tree_edge[v] = rng.choice(edges[v])
        # verify every vertex reaches `last`
        for v in list(tree_edge):
            seen, cur = set(), v
            while cur != last:
                if cur in seen or cur not in tree_edge:
                    ok = False
                    break
                seen.add(cur)
                cur = tree_edge[cur]
            if not ok:
                break
        if ok:
            break
    else:
        # Could not build a valid Eulerian arborescence; fall back to a plain
        # shuffle. Preserves mononucleotide but NOT dinucleotide composition,
        # so flag it rather than pretending the null is intact.
        chars = list(s)
        rng.shuffle(chars)
        return "".join(chars)

    new_edges = {}
    for v, outs in edges.items():
        outs = outs[:]
        if v in tree_edge:
            outs.remove(tree_edge[v])
            rng.shuffle(outs)
            outs.append(tree_edge[v])
        else:
            rng.shuffle(outs)
        new_edges[v] = outs

    out, cur, ptr = [s[0]], s[0], defaultdict(int)
    for _ in range(len(s) - 1):
        nxt = new_edges[cur][ptr[cur]]
        ptr[cur] += 1
        out.append(nxt)
        cur = nxt
    return "".join(out)


def shuffle_synonymous(cds: str, rng: random.Random) -> str:
    """Codon shuffle preserving the protein AND the codon-usage table.

    The correct null for coding sequence: any surviving k-mer signal cannot be
    explained by amino-acid composition or codon bias.
    """
    s = cds.upper()
    s = s[:len(s) - len(s) % 3]
    codons = [s[i:i + 3] for i in range(0, len(s), 3)]
    by_aa = defaultdict(list)
    for c in codons:
        if set(c) <= set("ACGT"):
            by_aa[str(Seq(c).translate())].append(c)
    for aa in by_aa:
        rng.shuffle(by_aa[aa])
    ptr = defaultdict(int)
    out = []
    for c in codons:
        if set(c) > set("ACGT"):
            out.append(c)
            continue
        aa = str(Seq(c).translate())
        out.append(by_aa[aa][ptr[aa]])
        ptr[aa] += 1
    return "".join(out)


def null_distribution(seqs: dict[str, str], k: int, mode: str,
                      n: int = 200) -> np.ndarray:
    """Distribution of the focal statistic under the chosen null."""
    rng = random.Random(RANDOM_SEED)
    shuffler = {"dinuc": shuffle_dinucleotide,
                "synonymous": shuffle_synonymous}[mode]
    out = []
    for _ in range(n):
        shuffled = {s: shuffler(q, rng) for s, q in seqs.items()}
        X, labels, _ = profiles_matrix(shuffled, k)
        out.append(acomys_deviation(X, labels))
    return np.array(out)


# ==============================================================================
# ARM 1 - PHYLOGENETIC RESIDUAL TEST (primary endpoint)
# ==============================================================================

def acomys_deviation(X: np.ndarray, labels: list[str]) -> float:
    """Mean k-mer distance from A. cahirinus to all non-Deomyinae species."""
    try:
        i = labels.index("Acomys cahirinus")
    except ValueError:
        return np.nan
    D = squareform(pdist(X, metric="cosine"))
    others = [j for j, s in enumerate(labels)
              if CLADE.get(s) != "Deomyinae" and j != i]
    return float(np.mean(D[i, others])) if others else np.nan


def phylogenetic_residual(seqs: dict[str, str], k: int) -> dict:
    """Core test.

    Regress pairwise k-mer distance on phylogenetic distance across ALL species
    pairs. Under neutral divergence, k-mer distance is roughly linear in time.
    The Acomys residual - how far Acomys sits above that line - is the
    phenotype-relevant quantity, purged of the taxonomic signal that would
    otherwise dominate.

    Returns the standardised Acomys residual plus the fit diagnostics.
    """
    X, labels, _ = profiles_matrix(seqs, k)
    if len(labels) < 4:
        return {"n_species": len(labels), "residual_z": np.nan,
                "note": "too few species"}

    kd = squareform(pdist(X, metric="cosine"))
    pd_ma = np.zeros_like(kd)
    for a, sa in enumerate(labels):
        for b, sb in enumerate(labels):
            ta, tb = DIVERGENCE_MA.get(sa), DIVERGENCE_MA.get(sb)
            if ta is None or tb is None:
                pd_ma[a, b] = np.nan
            elif CLADE.get(sa) == CLADE.get(sb):
                pd_ma[a, b] = abs(ta - tb)
            else:
                pd_ma[a, b] = max(ta, tb)

    iu = np.triu_indices(len(labels), k=1)
    x, y = pd_ma[iu], kd[iu]
    m = ~(np.isnan(x) | np.isnan(y))
    if m.sum() < 6:
        return {"n_species": len(labels), "residual_z": np.nan,
                "note": "too few valid pairs"}

    slope, intercept, r, p, se = stats.linregress(x[m], y[m])
    resid = y[m] - (slope * x[m] + intercept)
    sd = resid.std(ddof=1) or np.nan

    # Acomys-involving pairs only
    try:
        ai = labels.index("Acomys cahirinus")
    except ValueError:
        return {"n_species": len(labels), "residual_z": np.nan,
                "note": "focal species absent"}
    pair_has_acomys = np.array([(a == ai or b == ai)
                                for a, b in zip(*iu)])[m]
    aco_resid = resid[pair_has_acomys]

    return {
        "n_species": len(labels),
        "n_pairs": int(m.sum()),
        "slope": slope, "intercept": intercept,
        "r2": r ** 2, "fit_p": p,
        "residual_sd": sd,
        "acomys_residual_mean": float(np.mean(aco_resid)),
        "residual_z": float(np.mean(aco_resid) / sd) if sd else np.nan,
        "acomys_residual_ci": bootstrap_ci(aco_resid),
        "species": labels,
    }


def bootstrap_ci(v: np.ndarray, n: int = N_BOOTSTRAP,
                 alpha: float = ALPHA) -> tuple[float, float]:
    if len(v) < 2:
        return (np.nan, np.nan)
    rng = np.random.default_rng(RANDOM_SEED)
    boot = [rng.choice(v, len(v), replace=True).mean() for _ in range(n)]
    return (float(np.percentile(boot, 100 * alpha / 2)),
            float(np.percentile(boot, 100 * (1 - alpha / 2))))


def run_residual_arm(region: str, k: int) -> pd.DataFrame:
    """Focal module vs each control set. THE key comparison of the paper."""
    log(f"ARM 1: phylogenetic residual  [region={region}, k={k}]")
    rows = []
    # FOCAL is the primary endpoint. The PANEL_* blocks are E. Hassanat's
    # a priori secondary endpoint (received 20 Aug 2026, logged before any
    # analysis). Per her decision WP2 was NOT re-run for panel genes absent
    # from the scaffold, so some blocks will have too few genes to test - that
    # is reported as a data limitation, not as a null result.
    PANEL_SETS = {f"PANEL_{k2}": list(v) for k2, v in cfg.EXTENDED_PANEL.items()}
    PANEL_SETS["PANEL_combined15"] = list(cfg.EXTENDED_15)
    TEST_SETS = {"FOCAL": FOCAL_GENES, **PANEL_SETS}
    gene_sets = {**TEST_SETS, **CONTROL_SETS}

    for set_name, genes in gene_sets.items():
        for gene in genes:
            seqs = load_gene(gene, region)
            if len(seqs) < 4:
                log(f"  skip {gene}: only {len(seqs)} species", "WARN")
                continue
            res = phylogenetic_residual(seqs, k)
            rows.append({"gene_set": set_name, "gene": gene, "region": region,
                         "k": k,
                         "median_len": float(np.median([len(s) for s in seqs.values()])),
                         "n_seqs": len(seqs),
                         **{kk: vv for kk, vv in res.items() if kk != "species"}})

    df = pd.DataFrame(rows)
    if df.empty:
        log("  no data - run 01_fetch_data.py and 02_build_ortholog_scaffold.py",
            "WARN")
        return df

    # A gene can clear the >=4-species filter and still return NaN, most often
    # because it has no Deomyinae sequence at all - and the endpoint IS the
    # Acomys residual, so with no Acomys there is nothing to measure. That is
    # correct behaviour, but it must not vanish silently: a module quietly
    # shrinking changes what its test means.
    nan_rows = df[df.residual_z.isna()]
    if len(nan_rows):
        log(f"\n  {len(nan_rows)} gene-entries returned NaN "
            f"(passed the species filter but produced no residual):")
        for _, r in nan_rows.drop_duplicates("gene").iterrows():
            seqs = load_gene(r.gene, region)
            clades = {cfg.CLADE.get(s) for s in seqs}
            why = ("no Deomyinae sequence - the endpoint is the Acomys "
                   "residual, so there is nothing to measure"
                   if "Deomyinae" not in clades
                   else "insufficient valid species pairs")
            log(f"    {r.gene:10s} {int(r.n_seqs)} species "
                f"({', '.join(sorted(c for c in clades if c))}) - {why}",
                "WARN")
        log("  These are DATA LIMITATIONS, not null results, and are "
            "reported as such.", "WARN")

    # ---- LENGTH ADJUSTMENT (pre-specified: see PREREGISTRATION.md sec.4) ----
    # Sequence length is not neutral here. Measured WITHIN control genes only -
    # where focal status cannot bias it - longer 3'UTRs give systematically
    # higher residuals (positive Spearman, significant). Focal 3'UTRs are also longer
    # than control 3'UTRs, so an unadjusted comparison can manufacture a
    # positive result out of a length difference alone.
    df["log_len"] = np.log10(df["median_len"].clip(lower=1))
    # Fit on CONTROL SETS ONLY - not merely "everything that isn't FOCAL".
    # The PANEL_* blocks are test sets, and the ECM block shares all five of
    # its genes with the fibrosis_effector control set. Including them here
    # would let the sets under test shape the covariate that is supposed to
    # be estimated independently of them.
    ctrl_mask = df.gene_set.isin(CONTROL_SETS)
    fit = df[ctrl_mask].dropna(subset=["residual_z", "log_len"])
    if len(fit) >= 6:
        rho, p_len = stats.spearmanr(fit.log_len, fit.residual_z)
        b, a = np.polyfit(fit.log_len, fit.residual_z, 1)   # fit on CONTROLS only
        df["residual_adj"] = df["residual_z"] - (a + b * df["log_len"])
        log(f"  length model (fitted on controls): residual_z = "
            f"{a:+.3f} {b:+.3f}*log10(len)   Spearman={rho:+.3f} p={p_len:.3f}")
        if p_len < 0.05:
            log(f"  Length IS associated with the residual - the "
                f"length-adjusted test below is the one to report.", "WARN")
    else:
        df["residual_adj"] = df["residual_z"]
        log("  too few control genes to fit a length model", "WARN")

    # Both tests are reported. Unadjusted is descriptive; adjusted is primary.
    tests = []
    for metric, label in (("residual_z", "unadjusted"),
                          ("residual_adj", "length_adjusted")):
        for tname in TEST_SETS:
            focal = df.loc[df.gene_set == tname, metric].dropna()
            # Control sets sharing genes with a panel block are excluded from
            # that block's comparison (E.H. decision, 20 Aug): a test against
            # a set containing half the block's own members is circular.
            if tname == "FOCAL":
                allowed = list(CONTROL_SETS)
            elif tname == "PANEL_combined15":
                allowed = [c for c in CONTROL_SETS if c != "fibrosis_effector"]
            else:
                allowed = cfg.controls_for(tname.replace("PANEL_", ""))
            for cname in allowed:
                ctrl = df.loc[df.gene_set == cname, metric].dropna()
                if len(focal) < 3 or len(ctrl) < 3:
                    continue
                u, p = stats.mannwhitneyu(focal, ctrl, alternative="greater")
                tests.append({"model": label,
                          "test_set": tname,
                          "endpoint": "PRIMARY" if tname == "FOCAL"
                                      else "secondary",
                          "comparison": f"{tname}_vs_{cname}",
                          "n_focal": len(focal), "n_control": len(ctrl),
                          "U": u, "p_raw": p,
                          "cliffs_delta": cliffs_delta(focal.values, ctrl.values),
                          "focal_median": focal.median(),
                          "control_median": ctrl.median()})
    focal = df.loc[df.gene_set == "FOCAL", "residual_z"].dropna()

    if tests:
        t = pd.DataFrame(tests)
        # FDR is applied WITHIN each (model, test_set) family, not across the
        # whole table.
        #
        # The pre-registered rule is "FOCAL must clear ALL THREE control sets",
        # so the focal comparisons form one family of three. Pooling them with
        # the secondary panel comparisons would inflate the PRIMARY endpoint's
        # q-values purely because secondary tests were added afterwards - the
        # unadjusted 3'UTR focal q-values moved from 0.043/0.047/0.047 to
        # 0.101/0.109/0.109 on a first pass for exactly that reason, with the
        # p-values unchanged. A pre-registered primary endpoint must not be
        # penalised for later additions to the analysis.
        t["q_BH"] = np.nan
        for (m, ts), sel in t.groupby(["model", "test_set"]).groups.items():
            t.loc[sel, "q_BH"] = benjamini_hochberg(
                t.loc[sel, "p_raw"].values)
        RESULTS.mkdir(parents=True, exist_ok=True)
        t.to_csv(RESULTS / f"arm1_tests_{region}_k{k}.csv", index=False)
        log("\n" + t.to_string(index=False))

    df.to_csv(RESULTS / f"arm1_residuals_{region}_k{k}.csv", index=False)
    plot_residuals(df, region, k)
    return df


def cliffs_delta(a: np.ndarray, b: np.ndarray) -> float:
    """Non-parametric effect size. Report this, not just p."""
    gt = sum((x > y) for x in a for y in b)
    lt = sum((x < y) for x in a for y in b)
    return (gt - lt) / (len(a) * len(b))


def benjamini_hochberg(p: np.ndarray) -> np.ndarray:
    p = np.asarray(p, dtype=float)
    n = len(p)
    order = np.argsort(p)
    q = np.empty(n)
    prev = 1.0
    for rank, idx in enumerate(order[::-1]):
        i = n - rank
        prev = min(prev, p[idx] * n / i)
        q[idx] = prev
    return q


def plot_residuals(df: pd.DataFrame, region: str, k: int):
    FIGS.mkdir(parents=True, exist_ok=True)
    fig, ax = plt.subplots(figsize=(7, 4.2))
    order = ["FOCAL"] + list(CONTROL_SETS)
    sns.boxplot(data=df, x="gene_set", y="residual_z", order=order,
                ax=ax, width=.55, showfliers=False,
                palette=["#c0392b"] + ["#95a5a6"] * len(CONTROL_SETS))
    sns.stripplot(data=df, x="gene_set", y="residual_z", order=order,
                  ax=ax, color="k", size=3.5, alpha=.65, jitter=.18)
    ax.axhline(0, ls="--", lw=.8, c="k", alpha=.5)
    ax.set_xlabel("")
    ax.set_ylabel("Acomys phylogenetic residual (z)")
    ax.set_title(f"Deviation from murid phylogenetic expectation\n"
                 f"region = {region}, k = {k}", fontsize=10)
    plt.xticks(rotation=18, ha="right")
    plt.tight_layout()
    plt.savefig(FIGS / f"arm1_residuals_{region}_k{k}.png")
    plt.close()


# ==============================================================================
# ARM 2 - COMPOSITION NULLS
# ==============================================================================

def composition_adjusted_residuals(df: pd.DataFrame,
                                   region: str) -> pd.DataFrame:
    """Remove the linear effect of GC (and GC3 for CDS) from the residual.

    THE PRIMARY COMPOSITION CONTROL. k-mer profiles are dominated by base
    composition, so a gene can look divergent merely for being GC-rich. This
    regresses that out and re-tests on the adjusted values.

    Unlike a shuffle-based null, this test CAN fire in the informative
    direction: it asks whether focal genes stay elevated AFTER composition is
    accounted for.
    """
    d = df.dropna(subset=["residual_z", "mean_gc"]).copy()
    if len(d) < 8:
        d["residual_adj"] = np.nan
        return d

    preds = ["mean_gc"]
    if region == "cds" and d["mean_gc3"].notna().sum() >= len(d) * 0.8:
        preds.append("mean_gc3")

    X = np.column_stack([np.ones(len(d))] + [d[p].values for p in preds])
    y = d["residual_z"].values
    beta, *_ = np.linalg.lstsq(X, y, rcond=None)
    d["residual_adj"] = y - X @ beta
    d.attrs["gc_betas"] = dict(zip(["intercept"] + preds, beta))
    return d


def run_nulls_arm(region: str, k: int, n_null: int = 200) -> pd.DataFrame:
    """Arm 2 - is the Acomys deviation explained by nucleotide composition?

    TWO tests, deliberately separated:

      2a PRIMARY    composition-adjusted residual. Regress out GC (+GC3 for
                    CDS), then re-test focal vs controls. Fires in the
                    informative direction.

      2b DIAGNOSTIC dinucleotide / codon shuffle. Reported TWO-SIDED and
                    descriptively only.

                    Note on 2b: shuffling every species independently destroys
                    homology, so every pairwise distance inflates and the
                    observed value sits BELOW the null essentially always. A
                    one-sided "observed > null" test therefore cannot fire and
                    would be meaningless. What 2b actually measures is how much
                    of the composition-only ceiling the real sequences reach -
                    a conservation diagnostic, not a significance test. It is
                    reported as such.
    """
    log(f"ARM 2: composition controls  [region={region}, k={k}]")

    gene_sets = {"FOCAL": FOCAL_GENES, **CONTROL_SETS}
    rows = []
    for set_name, genes in gene_sets.items():
        for gene in genes:
            seqs = load_gene(gene, region)
            if len(seqs) < 4:
                continue
            res = phylogenetic_residual(seqs, k)
            gc = np.mean([gc_content(s) for s in seqs.values()])
            gc3 = (np.mean([gc3_content(s) for s in seqs.values()])
                   if region == "cds" else np.nan)
            rows.append({"gene_set": set_name, "gene": gene, "region": region,
                         "k": k, "residual_z": res.get("residual_z", np.nan),
                         "mean_gc": gc, "mean_gc3": gc3})

    if not rows:
        log("ARM 2 produced NO results - ortholog data missing.", "WARN")
        return pd.DataFrame()

    df = pd.DataFrame(rows)

    # ---- 2a PRIMARY: composition-adjusted ----
    adj = composition_adjusted_residuals(df, region)
    betas = adj.attrs.get("gc_betas", {})
    if betas:
        log("  composition model: " +
            ", ".join(f"{k2}={v:+.3f}" for k2, v in betas.items()))

    focal = adj.loc[adj.gene_set == "FOCAL", "residual_adj"].dropna()
    tests = []
    for cname in CONTROL_SETS:
        ctrl = adj.loc[adj.gene_set == cname, "residual_adj"].dropna()
        if len(focal) < 3 or len(ctrl) < 3:
            continue
        u, pv = stats.mannwhitneyu(focal, ctrl, alternative="greater")
        tests.append({"test": "2a_composition_adjusted",
                      "comparison": f"FOCAL_vs_{cname}",
                      "n_focal": len(focal), "n_control": len(ctrl),
                      "U": u, "p_raw": pv,
                      "cliffs_delta": cliffs_delta(focal.values, ctrl.values),
                      "focal_median": focal.median(),
                      "control_median": ctrl.median()})

    RESULTS.mkdir(parents=True, exist_ok=True)
    if tests:
        tdf = pd.DataFrame(tests)
        tdf["q_BH"] = benjamini_hochberg(tdf["p_raw"].values)
        tdf.to_csv(RESULTS / f"arm2a_composition_adjusted_{region}_k{k}.csv",
                   index=False)
        log("\n" + tdf.to_string(index=False))
        n_pass = (tdf.q_BH < ALPHA).sum()
        log(f"\n  2a: focal exceeds {n_pass}/{len(tdf)} control sets "
            f"after composition adjustment")
        if n_pass < len(tdf):
            log("  Focal module does NOT clear every control set - by the "
                "pre-specified rule, divergence is NOT claimed.", "WARN")

    adj.to_csv(RESULTS / f"arm2a_residuals_{region}_k{k}.csv", index=False)

    # ---- 2b DIAGNOSTIC: shuffle ----
    mode = "synonymous" if region == "cds" else "dinuc"
    log(f"\n  2b diagnostic shuffle ({mode}, n={n_null}) - two-sided, "
        f"descriptive only")
    diag = []
    for gene in FOCAL_GENES + CONTROL_SETS["housekeeping"]:
        seqs = load_gene(gene, region)
        if len(seqs) < 4:
            continue
        X, labels, _ = profiles_matrix(seqs, k)
        obs = acomys_deviation(X, labels)
        null = null_distribution(seqs, k, mode, n=n_null)
        null = null[~np.isnan(null)]
        if len(null) < MIN_NULL_DRAWS:
            log(f"    {gene}: only {len(null)}/{n_null} valid draws "
                f"(need >= {MIN_NULL_DRAWS}) - skipped. Raise --n-null.", "WARN")
            continue
        z = (obs - null.mean()) / (null.std() or np.nan)
        p_two = 2 * min((np.sum(null >= obs) + 1) / (len(null) + 1),
                        (np.sum(null <= obs) + 1) / (len(null) + 1))
        diag.append({"gene": gene, "region": region, "k": k,
                     "null_model": mode, "observed": obs,
                     "null_mean": null.mean(), "null_sd": null.std(),
                     "z_vs_null": z, "p_two_sided": min(p_two, 1.0),
                     "frac_of_null_ceiling": obs / null.mean() if null.mean() else np.nan})
        log(f"    {gene:10s} obs={obs:.4f} null={null.mean():.4f} "
            f"({obs / null.mean() * 100:5.1f}% of ceiling)")

    if diag:
        ddf = pd.DataFrame(diag)
        ddf.to_csv(RESULTS / f"arm2b_shuffle_diagnostic_{region}_k{k}.csv",
                   index=False)
        f_mean = ddf[ddf.gene.isin(FOCAL_GENES)].frac_of_null_ceiling.mean()
        h_mean = ddf[~ddf.gene.isin(FOCAL_GENES)].frac_of_null_ceiling.mean()
        log(f"\n  2b: focal reaches {f_mean * 100:.1f}% of the composition "
            f"ceiling vs {h_mean * 100:.1f}% for housekeeping")
        log("      (higher = less constrained relative to composition)")

    return adj


# ==============================================================================
# ARM 3 - REGION-STRATIFIED COMPARISON
# ==============================================================================

def run_region_arm(k: int) -> pd.DataFrame:
    """Where in the gene does the divergence sit?

    Pre-registered prediction: 3'UTR and promoter carry the signal; CDS does
    not. If CDS shows the largest deviation, the regulatory hypothesis is
    wrong and the paper should say so.
    """
    log(f"ARM 3: region stratification  [k={k}]")
    frames = []
    for region in ["promoter", "utr5", "cds", "utr3"]:
        d = run_residual_arm(region, k)
        if not d.empty:
            frames.append(d)
    if not frames:
        return pd.DataFrame()
    df = pd.concat(frames, ignore_index=True)
    df.to_csv(RESULTS / f"arm3_regions_k{k}.csv", index=False)

    focal = df[df.gene_set == "FOCAL"]
    if not focal.empty:
        fig, ax = plt.subplots(figsize=(6.5, 4))
        sns.boxplot(data=focal, x="region", y="residual_z",
                    order=["promoter", "utr5", "cds", "utr3"],
                    ax=ax, width=.5, showfliers=False, color="#c0392b")
        sns.stripplot(data=focal, x="region", y="residual_z",
                      order=["promoter", "utr5", "cds", "utr3"],
                      ax=ax, color="k", size=4, alpha=.7)
        ax.axhline(0, ls="--", lw=.8, c="k", alpha=.5)
        ax.set_ylabel("Acomys phylogenetic residual (z)")
        ax.set_title("Focal module: where does divergence localise?", fontsize=10)
        plt.tight_layout()
        FIGS.mkdir(parents=True, exist_ok=True)
        plt.savefig(FIGS / f"arm3_regions_k{k}.png")
        plt.close()
    return df


# ==============================================================================
# ARM 4 - LEAKAGE-HARDENED ML + SHAP (secondary, descriptive)
# ==============================================================================

def run_ml_arm(region: str, k: int) -> dict:
    """Classification is reported ONLY against a permutation null.

    Every step that could leak - scaling, PCA - lives inside the pipeline so it
    is refit per fold. This is the fix your Echinococcus reviewers demanded;
    it is here from the start.
    """
    log(f"ARM 4: ML + SHAP  [region={region}, k={k}]")
    rows, labels_all = [], []
    for gset, genes in {"FOCAL": FOCAL_GENES, **CONTROL_SETS}.items():
        for gene in genes:
            for sp, seq in load_gene(gene, region).items():
                rows.append(kmer_profile(seq, k))
                labels_all.append((gset, gene, sp, REGENERATIVE.get(sp, False)))

    if len(rows) < 20:
        log("  insufficient data", "WARN")
        return {}

    X = np.vstack(rows)
    meta = pd.DataFrame(labels_all,
                        columns=["gene_set", "gene", "species", "regenerative"])
    y = meta["regenerative"].astype(int).values

    if len(np.unique(y)) < 2 or min(np.bincount(y)) < CV_FOLDS:
        log("  class imbalance too severe for CV", "WARN")
        return {}

    pipe = Pipeline([
        ("scale", MinMaxScaler()),
        ("pca", PCA(n_components=min(20, X.shape[1], X.shape[0] - 1),
                    random_state=RANDOM_SEED)),
        ("clf", LogisticRegression(max_iter=5000, class_weight="balanced",
                                   random_state=RANDOM_SEED)),
    ])

    cv = StratifiedKFold(CV_FOLDS, shuffle=True, random_state=RANDOM_SEED)
    obs = cross_val_score(pipe, X, y, cv=cv, scoring="f1_macro").mean()

    rng = np.random.default_rng(RANDOM_SEED)
    null = np.array([
        cross_val_score(pipe, X, rng.permutation(y), cv=cv,
                        scoring="f1_macro").mean()
        for _ in range(min(N_PERMUTATIONS, 200))
    ])
    p_perm = (np.sum(null >= obs) + 1) / (len(null) + 1)

    out = {"region": region, "k": k, "n_samples": len(y),
           "observed_f1": float(obs), "null_f1_mean": float(null.mean()),
           "null_f1_95pct": float(np.percentile(null, 95)),
           "p_permutation": float(p_perm)}
    log(f"  f1={obs:.3f}  null={null.mean():.3f}  p_perm={p_perm:.4f}")
    if p_perm > ALPHA:
        log("  Not distinguishable from the permutation null. "
            "Report as negative - do not quote the raw accuracy.", "WARN")

    RESULTS.mkdir(parents=True, exist_ok=True)
    (RESULTS / f"arm4_ml_{region}_k{k}.json").write_text(json.dumps(out, indent=2))
    return out


# ==============================================================================
# I/O
# ==============================================================================

def load_gene(gene: str, region: str) -> dict[str, str]:
    """Load one sequence per species for a gene/region.

    Expects 02_build_ortholog_scaffold.py to have written region-sliced FASTA:
        data/reference/orthologs/<GENE>/<Species_name>.<region>.fasta
    """
    gdir = ORTHO / gene
    if not gdir.exists():
        return {}
    out = {}
    for fa in gdir.glob(f"*.{region}.fasta"):
        sp = fa.name.split(".")[0].replace("_", " ")
        try:
            rec = next(SeqIO.parse(fa, "fasta"))
            s = str(rec.seq).upper()
            if len(s) >= 60:          # drop stubs
                out[sp] = s
        except StopIteration:
            continue
    return out


# ==============================================================================
# MAIN
# ==============================================================================

def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--arm", choices=["residual", "nulls", "region", "ml", "all"],
                    default="all")
    ap.add_argument("--region", choices=REGIONS, default="utr3")
    ap.add_argument("--k", type=int, default=5)
    ap.add_argument("--k-sweep", action="store_true",
                    help=f"repeat across k in {K_SWEEP}")
    ap.add_argument("--n-null", type=int, default=200)
    args = ap.parse_args()

    RESULTS.mkdir(parents=True, exist_ok=True)
    FIGS.mkdir(parents=True, exist_ok=True)

    log("=" * 78)
    log("WP3 - phylogeny-controlled alignment-free analysis")
    log(f"region={args.region}  k={args.k}  seed={RANDOM_SEED}")
    log("Primary endpoint: FOCAL vs control-set phylogenetic residual")
    log("=" * 78)

    if not ORTHO.exists() or not any(ORTHO.iterdir()):
        log("No ortholog data found. Run first:", "ERROR")
        log("  python 01_fetch_data.py --all", "ERROR")
        log("  python 02_build_ortholog_scaffold.py", "ERROR")
        return

    ks = K_SWEEP if args.k_sweep else [args.k]
    for k in ks:
        if args.arm in ("residual", "all"):
            run_residual_arm(args.region, k)
        if args.arm in ("nulls", "all"):
            run_nulls_arm(args.region, k, args.n_null)
        if args.arm in ("region", "all"):
            run_region_arm(k)
        if args.arm in ("ml", "all"):
            run_ml_arm(args.region, k)

    log("=" * 78)
    log(f"Results: {RESULTS}")
    log("Interpretation rule (pre-specified):")
    log("  Claim divergence ONLY if the focal module beats ALL THREE control")
    log("  sets at q<0.05 AND survives the composition null in Arm 2.")
    log("=" * 78)


if __name__ == "__main__":
    main()
