#!/usr/bin/env python3
"""
================================================================================
ALIGNMENT INSPECTION - run between 05 --align and 05 --run
================================================================================

--align reports how many codons survive codeml's cleandata=1 per GENE. That is
the number that decides whether a gene is analysable. It is not the number that
tells you WHY a gene lost codons, and the why matters more.

The failure mode this script exists to catch:

    If ONE species is systematically gappy, it deletes columns from EVERY gene
    it appears in. Column deletion is all-or-nothing across species, so a
    single bad sequence silently sets the ceiling for the whole analysis.

    That is bad in general. It is far worse if the bad species sits in the
    FOREGROUND clade, because branch-site model A estimates a separate omega
    on exactly that branch. Poor sequence on the foreground branch is the
    textbook route to a spurious positive selection result - alignment error
    is read by codeml as nonsynonymous substitution.

A. cahirinus is the specific worry here: its annotation is CAT-projected from
GRCm39 (see the acomys_annotation_circularity design gap), and it is a
foreground species.

Usage:
    python code\\inspect_alignments.py
    python code\\inspect_alignments.py --gene CX3CR1     # dump one alignment
================================================================================
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import config as cfg

try:
    from Bio import SeqIO
except ImportError:
    print("biopython required:  pip install biopython")
    sys.exit(1)

log = cfg.log
ALN = cfg.RESULTS / "wp4_selection" / "alignments"
OUT = cfg.RESULTS / "wp4_selection"
FOREGROUND_CLADE = "Deomyinae"
MIN_CLEAN_CODONS = 100


def gene_sets() -> dict[str, str]:
    """gene -> set name."""
    d = {g: "FOCAL" for g in cfg.FOCAL_GENES}
    for k, v in cfg.CONTROL_SETS.items():
        for g in (v or []):
            d.setdefault(g, k)
    return d


def read_aln(p: Path) -> dict[str, str]:
    return {r.id.replace("_", " "): str(r.seq).upper()
            for r in SeqIO.parse(p, "fasta")}


def codon_ok(c: str) -> bool:
    return len(c) == 3 and set(c) <= set("ACGT")


def analyse(seqs: dict[str, str]) -> tuple[int, int, dict[str, float]]:
    """Return (n_codons, n_clean, per-species bad-codon fraction)."""
    if not seqs:
        return 0, 0, {}
    L = min(len(s) for s in seqs.values()) // 3
    bad = {sp: 0 for sp in seqs}
    clean = 0
    for i in range(L):
        cols = {sp: s[i * 3:i * 3 + 3] for sp, s in seqs.items()}
        ok = True
        for sp, c in cols.items():
            if not codon_ok(c):
                bad[sp] += 1
                ok = False
        clean += ok
    return L, clean, {sp: (b / L if L else 0.0) for sp, b in bad.items()}


def main(dump: str | None):
    if not ALN.exists():
        log(f"not found: {ALN}  - run 05_selection_analysis.py --align", "ERROR")
        return

    files = sorted(ALN.glob("*.codon.aln"))
    if not files:
        log(f"no *.codon.aln in {ALN}", "ERROR")
        return

    sets = gene_sets()
    rows = []
    bad_by_sp: dict[str, list[float]] = {}
    seen_by_sp: dict[str, int] = {}

    for f in files:
        gene = f.name.replace(".codon.aln", "")
        seqs = read_aln(f)
        L, clean, bad = analyse(seqs)
        rows.append({
            "gene": gene,
            "set": sets.get(gene, "?"),
            "n_sp": len(seqs),
            "codons": L,
            "clean": clean,
            "pct_kept": round(100 * clean / L, 1) if L else 0.0,
            "usable": clean >= MIN_CLEAN_CODONS,
        })
        for sp, frac in bad.items():
            bad_by_sp.setdefault(sp, []).append(frac)
            seen_by_sp[sp] = seen_by_sp.get(sp, 0) + 1

    if dump:
        f = ALN / f"{dump}.codon.aln"
        if not f.exists():
            log(f"no alignment for {dump}", "ERROR")
            return
        seqs = read_aln(f)
        L, clean, bad = analyse(seqs)
        log(f"\n{dump}: {len(seqs)} species, {L} codons, {clean} clean")
        log("\nper-species codons that are gapped/ambiguous:")
        for sp, fr in sorted(bad.items(), key=lambda x: -x[1]):
            mark = "  <-- FOREGROUND" if cfg.CLADE.get(sp) == FOREGROUND_CLADE else ""
            log(f"  {sp:26s} {fr*100:5.1f}%{mark}")
        log("\nfirst 40 codons, translated:")
        from Bio.Seq import Seq
        for sp, s in seqs.items():
            pep = str(Seq(s[:120]).translate())
            log(f"  {sp:26s} {pep}")
        return

    df = pd.DataFrame(rows).sort_values(["set", "clean"])
    log("=" * 74)
    log("PER-GENE ALIGNMENT QUALITY")
    log("=" * 74)
    log("\n" + df.to_string(index=False))
    df.to_csv(OUT / "alignment_qc.csv", index=False)

    log("\n" + "=" * 74)
    log("SURVIVAL BY GENE SET  (this is what the Fisher test compares)")
    log("=" * 74)
    g = df.groupby("set").agg(n=("gene", "size"), usable=("usable", "sum"),
                              median_clean=("clean", "median"),
                              median_pct=("pct_kept", "median"))
    g["pct_usable"] = (100 * g.usable / g.n).round(1)
    log("\n" + g.to_string())

    if len(g) > 1:
        spread = g.pct_usable.max() - g.pct_usable.min()
        if spread > 25:
            log(f"\nWARNING: usable-gene rates differ by {spread:.0f} points "
                "across sets. A focal-vs-control difference in sequence "
                "QUALITY would masquerade as one in selection.", "WARN")
        else:
            log(f"\nOK: usable-gene rates within {spread:.0f} points across "
                "sets - the comparison is like-for-like.")

    log("\n" + "=" * 74)
    log("PER-SPECIES COLUMN LOSS  (the ceiling-setter)")
    log("=" * 74)
    log("\nMean % of codons that are gapped/ambiguous, averaged over genes.")
    log("A high value here costs columns in EVERY gene that species appears in.\n")
    sp_rows = []
    for sp, fracs in bad_by_sp.items():
        sp_rows.append({"species": sp,
                        "clade": cfg.CLADE.get(sp, "?"),
                        "n_genes": seen_by_sp[sp],
                        "mean_bad_pct": round(100 * float(np.mean(fracs)), 1),
                        "max_bad_pct": round(100 * float(np.max(fracs)), 1)})
    sdf = pd.DataFrame(sp_rows).sort_values("mean_bad_pct", ascending=False)
    sdf["foreground"] = sdf.clade == FOREGROUND_CLADE
    log(sdf.to_string(index=False))
    sdf.to_csv(OUT / "alignment_qc_by_species.csv", index=False)

    worst = sdf.iloc[0]
    if worst.foreground:
        log(f"\nWARNING: the gappiest species is {worst.species}, which is in "
            f"the FOREGROUND clade ({worst.mean_bad_pct}% mean loss).", "WARN")
        log("Branch-site model A fits a separate omega on precisely this "
            "branch. Alignment error there is read as nonsynonymous change, "
            "which inflates foreground omega and manufactures positive "
            "selection. Treat any positive result as suspect until this is "
            "resolved - the leave-one-out check below is the minimum.", "WARN")
        log(f"\n  Sensitivity check:  drop {worst.species} and re-run. If the "
            "result survives without it, it is not an artefact of that "
            "sequence. If it does not, the result was that sequence.", "WARN")
    else:
        fg = sdf[sdf.foreground]
        if len(fg):
            log(f"\nOK: no foreground species is the worst offender. "
                f"Foreground mean loss {fg.mean_bad_pct.min():.1f}"
                f"-{fg.mean_bad_pct.max():.1f}%, vs "
                f"{worst.mean_bad_pct:.1f}% for {worst.species}.")

    log(f"\nwrote {OUT / 'alignment_qc.csv'}")
    log(f"      {OUT / 'alignment_qc_by_species.csv'}")
    log("\nNext: eyeball 2-3 alignments, focal ones first.")
    log("  python code\\inspect_alignments.py --gene CX3CR1")


def foreground_audit():
    """Why is each foreground species present or absent, gene by gene?

    Branch-site model A is only as good as the foreground branch. If foreground
    species are missing from some genes and not others - and especially if that
    pattern differs between focal and control sets - the test is comparing
    trees, not selection.
    """
    from Bio.Seq import Seq
    sets = gene_sets()
    genes = sorted(sets)
    fg_sp = [s for s, c in cfg.CLADE.items() if c == FOREGROUND_CLADE]

    def drop_reason(nt: str) -> str:
        s = "".join(c if c in "ACGTacgt" else "N" for c in nt).upper()
        if len(s) < 60:
            return "too short"
        t = len(s) % 3
        if t:
            s = s[:-t]
        pep = str(Seq(s).translate())
        if pep.endswith("*"):
            pep = pep[:-1]
        n = pep.count("*")
        if n:
            return f"frame broken ({n} internal stops)"
        if pep.count("X") / max(len(pep), 1) > 0.10:
            return "ambiguous"
        return "PASS"

    log("=" * 74)
    log("FOREGROUND AUDIT - why each Deomyinae species is in or out")
    log("=" * 74)

    rows = []
    for sp in fg_sp:
        key = sp.replace(" ", "_")
        for g in genes:
            f = cfg.ORTHO / g / f"{key}.cds.fasta"
            if not f.exists():
                r = "no CDS in scaffold"
            else:
                try:
                    seq = str(next(SeqIO.parse(f, "fasta")).seq)
                    r = drop_reason(seq)
                except StopIteration:
                    r = "empty file"
            aln = ALN / f"{g}.codon.aln"
            in_aln = aln.exists() and sp in read_aln(aln)
            rows.append({"species": sp, "gene": g, "set": sets[g],
                         "scaffold": r, "in_alignment": in_aln})

    df = pd.DataFrame(rows)
    df.to_csv(OUT / "foreground_audit.csv", index=False)

    log("\nScaffold status by species:")
    simple = df.scaffold.str.replace(r" \(.*\)", "", regex=True)
    log("\n" + pd.crosstab(df.species, simple).to_string())

    log("\nIn the ALIGNMENT, by species and gene set:")
    log("\n" + pd.crosstab([df.species, df.set], df.in_alignment).to_string())

    log("\nForeground tips per gene, by gene set:")
    per = (df[df.in_alignment].groupby(["set", "gene"]).size()
           .rename("n_fg").reset_index())
    allg = pd.DataFrame({"gene": genes, "set": [sets[g] for g in genes]})
    per = allg.merge(per, on=["gene", "set"], how="left").fillna({"n_fg": 0})
    log("\n" + pd.crosstab(per.set, per.n_fg.astype(int)).to_string())
    means = per.groupby("set").n_fg.mean()
    log("\nMean foreground tips per gene:")
    log("\n" + means.round(2).to_string())

    if "FOCAL" in means.index:
        ctrl = means.drop("FOCAL")
        gap = means["FOCAL"] - ctrl.max()
        if gap > 0.3:
            log(f"\nSTOP: the focal set averages {means['FOCAL']:.2f} "
                f"foreground tips per gene against {ctrl.max():.2f} for the "
                "best control set.", "ERROR")
            log("Branch-site power scales with foreground branch length, so "
                "focal genes are more likely to reach significance from tree "
                "structure alone - independent of any selection.", "ERROR")
            log("With 3 Acomys the foreground is an ANCESTRAL CLADE branch; "
                "with 1 it is a TERMINAL branch. Those are different "
                "hypotheses and must not be pooled in one Fisher test.",
                "ERROR")
            log("Fix: restrict every gene to one fixed species set so the "
                "topology and the foreground are identical throughout.",
                "ERROR")
        else:
            log(f"\nOK: focal and control sets carry comparable foreground "
                f"representation (gap {gap:.2f} tips).")

    log(f"\nwrote {OUT / 'foreground_audit.csv'}")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--gene", help="dump one alignment in detail")
    ap.add_argument("--foreground", action="store_true",
                    help="audit foreground species coverage (run this)")
    a = ap.parse_args()
    if a.foreground:
        foreground_audit()
    else:
        main(a.gene)
