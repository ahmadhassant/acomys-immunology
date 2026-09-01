#!/usr/bin/env python3
"""
================================================================================
WP4 - H2: BRANCH-SPECIFIC SELECTION ON THE ACOMYS LINEAGE
PAML codeml, branch-site model A
================================================================================

WHY THIS ARM IS WORTH RUNNING WHEN THREE OTHERS FAILED
-------------------------------------------------------
The spleen arm had no data. H1 was a sequence-LENGTH artefact. H4 was a FACS
SORT artefact. This arm inherits neither:

  * CDS recovery was 1.00-1.03x of query across all five focal genes - the
    cleanest data in the project (see WP2 ortholog_table.csv).
  * Coding sequence length is constrained by the reading frame, so the length
    confound that sank H1 cannot operate the same way here.
  * This is genomic sequence, not sorted cells, so no gating artefact exists.
  * dN/dS is a RATIO of two rates measured on the same alignment, so it is
    internally controlled against alignment quality and composition in a way
    k-mer distance is not.

It is also the only arm testing sequence EVOLUTION rather than sequence
COMPOSITION. H1 showed composition carries no signal; that says nothing about
whether selection acted.

THE TEST
--------
codeml branch-site model A, Acomys branch as foreground:

    Null       : model=2 NSsites=2 fix_omega=1 omega=1
    Alternative: model=2 NSsites=2 fix_omega=0
    LRT        : 2*(lnL_alt - lnL_null) ~ chi2, df=1
                 (a 50:50 mixture; using df=1 is the standard CONSERVATIVE
                  choice recommended by Yang & dos Reis 2011)

Applied IDENTICALLY to the focal module and all three control sets, with
BH-FDR across the whole panel. The comparison of interest is not "is any
focal gene under selection" - with 35 genes something usually is - but
"is the focal module enriched for branch-specific selection RELATIVE to
matched controls".

PRE-SPECIFIED DECISION RULE (mirrors H1's, fixed before running)
    Claim branch-specific selection ONLY IF the focal module shows a higher
    proportion of genes with significant LRT (q<0.05) than ALL THREE control
    sets, tested by Fisher's exact test.

Usage:
    python code\\05_selection_analysis.py --check
    python code\\05_selection_analysis.py --align        # codon-aware MSA
    python code\\05_selection_analysis.py --run          # codeml, slow
    python code\\05_selection_analysis.py --summarise

Requires: mafft, codeml (PAML). biopython.
================================================================================
"""

from __future__ import annotations

import argparse
import os
import re
import shutil
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from scipy import stats

sys.path.insert(0, str(Path(__file__).resolve().parent))
import config as cfg
import platform_compat as pc

cfg.set_all_seeds()
log = cfg.log

OUT = cfg.RESULTS / "wp4_selection"
ALN = OUT / "alignments"
RUNS = OUT / "codeml"

# Foreground = the Acomys (Deomyinae) lineage.
FOREGROUND_CLADE = "Deomyinae"
MIN_SPECIES_FOR_CODEML = 6       # below this the LRT is not worth running
MIN_CDS_LENGTH = 150

# ==============================================================================
# THE SPECIES PANEL IS FIXED ACROSS ALL GENES  (added 19 Aug 2026)
# ==============================================================================
# Letting each gene use whatever species happened to be recovered was a
# confound, not a convenience. Measured on the first alignment build:
#
#   A. russatus    34/35 genes      A. cahirinus   7/35      A. dimidiatus 7/35
#
#   mean foreground tips per gene:
#     FOCAL 2.2 | housekeeping 1.6 | immune 1.1 | fibrosis 1.0
#
# Two independent problems followed. Branch-site power scales with foreground
# branch length, so the focal set was more likely to reach significance from
# tree structure alone. And with 3 Acomys the '#1' mark sits on an ANCESTRAL
# CLADE branch while with 1 it sits on a TERMINAL branch - different
# hypotheses, which the Fisher test would have pooled.
#
# WHY THE OTHER TWO ACOMYS ARE EXCLUDED, and why re-running WP2 cannot fix it:
# A. cahirinus and A. dimidiatus are CAT-projected from GRCm39. Their CDS is
# frame-broken in 27/35 genes (internal stop codons). A. russatus, which is
# independently RefSeq-annotated (GCF_903995435.1), is clean in 34/35 with
# zero frame breaks. Testing all six reading frames and both strands recovered
# only 3 of the 27 - the remaining 24 carry SCATTERED INTERNAL INDELS, so no
# reframing rescues them. The projected coordinates do not respect the actual
# Acomys sequence. This is the acomys_annotation_circularity design gap
# manifesting as a hard blocker.
#
# WHAT THE SINGLE-SPECIES FOREGROUND ACTUALLY TESTS. Because A. russatus is the
# only Deomyinae in the tree, its terminal branch spans the whole Deomyinae
# lineage back to the Gerbillinae split (~18-20 Ma). Selection anywhere on that
# lineage - including the ancestral branch where regeneration presumably arose
# - falls on this branch. What is LOST is the ability to localise selection
# within the lineage: ancestral Deomyinae and A. russatus-specific events are
# not separable. State this in the manuscript; do not claim the ancestral
# branch specifically.
#
# A side benefit worth stating: the foreground is now the one Deomyinae genome
# NOT projected from mouse, so this arm is free of the annotation circularity
# that constrains the rest of the project.
H2_PANEL = [
    "Acomys russatus",          # Deomyinae  - FOREGROUND, RefSeq-annotated
    "Meriones unguiculatus",    # Gerbillinae - sister clade
    "Psammomys obesus",         # Gerbillinae
    "Apodemus sylvaticus",      # Murinae
    "Mastomys coucha",          # Murinae
    "Mus musculus",             # Murinae
    "Rattus norvegicus",        # Murinae
    "Mesocricetus auratus",     # Cricetidae - root
    "Peromyscus maniculatus",   # Cricetidae
]
# Genes missing any panel member are dropped, so every LRT uses one topology.
# 33/35 genes retain the complete panel (ACTB lacks M. coucha; CD3E lacks
# A. russatus and so has no foreground at all).
USE_FIXED_PANEL = True


# ==============================================================================
# TREE
# ==============================================================================
# Newick with the foreground branch marked '#1' as codeml requires.
# Topology from Steppan et al. 2004 / Alhajeri et al. 2015: Deomyinae and
# Gerbillinae are sisters, Murinae is outgroup to both, Cricetidae roots.

# Within-clade branching order, outermost (earliest-diverging) FIRST.
# A flat clade like (A,B,C,D) is a polytomy, and PAML's compiled MAXNSONS
# limit rejects it outright: "error: too many daughter nodes, raise MAXNSONS".
# Beyond that mechanical constraint, a polytomy is a claim of simultaneous
# divergence that nobody believes; resolving it is the honest default.
#
# Murinae order follows the tribal arrangement in Steppan & Schenk (2017) and
# Aghova et al. (2018): Rattini earliest, then Apodemyini, then Murini sister
# to Arvicanthini + Praomyini. The deeper Murinae nodes are NOT firmly
# resolved in the literature and this should be treated as one defensible
# topology rather than the settled answer.
#
# Why that uncertainty is tolerable here: the SAME topology is applied to
# focal and control genes alike, so topology error is non-differential between
# the sets being compared. It cannot generate a focal-vs-control difference.
# It could still bias absolute branch lengths, so an alternative-topology
# sensitivity run is worth reporting if H2 comes out positive.
CLADE_ORDER = {
    "Murinae": ["Rattus norvegicus", "Apodemus sylvaticus", "Mus musculus",
                "Mastomys coucha", "Grammomys surdaster"],
    # A. cahirinus and A. dimidiatus are near-conspecific; russatus is the
    # earlier split within the genus.
    "Deomyinae": ["Acomys russatus", "Acomys cahirinus", "Acomys dimidiatus"],
    "Gerbillinae": ["Meriones unguiculatus", "Psammomys obesus"],
    "Cricetidae": ["Mesocricetus auratus", "Peromyscus maniculatus"],
}


def node_degrees(newick: str) -> list[int]:
    """Children per internal node, root last. Guards the MAXNSONS limit."""
    stack: list[int] = []
    out: list[int] = []
    for ch in newick.strip().rstrip(";"):
        if ch == "(":
            stack.append(0)
        elif ch == ")":
            out.append(stack.pop() + 1)
        elif ch == "," and stack:
            stack[-1] += 1
    return out


def _comb(names: list[str]) -> str:
    """Nest a list into a strictly bifurcating comb: (a,(b,(c,d)))."""
    if len(names) == 1:
        return names[0]
    if len(names) == 2:
        return f"({names[0]},{names[1]})"
    return f"({names[0]},{_comb(names[1:])})"


def build_tree(species: list[str]) -> str:
    """Newick for the species present, with the Deomyinae clade as foreground.

    Strictly bifurcating within clades, trifurcating at the root (= unrooted).
    """
    by_clade: dict[str, list[str]] = {}
    for s in species:
        by_clade.setdefault(cfg.CLADE.get(s, "Murinae"), []).append(s)

    def grp(clade_names, clade=None):
        if not clade_names:
            return None
        order = CLADE_ORDER.get(clade, [])
        # Known species in their branching order, then any unknown ones
        # appended deterministically so the tree is reproducible.
        ordered = ([s for s in order if s in clade_names]
                   + sorted(s for s in clade_names if s not in order))
        return _comb([s.replace(" ", "_") for s in ordered])

    deo = grp(by_clade.get("Deomyinae", []), "Deomyinae")
    ger = grp(by_clade.get("Gerbillinae", []), "Gerbillinae")
    mur = grp(by_clade.get("Murinae", []), "Murinae")
    cri = grp(by_clade.get("Cricetidae", []), "Cricetidae")

    if deo is None:
        return ""
    deo = deo + " #1"                       # foreground mark

    # codeml requires an UNROOTED tree. A bifurcating root adds a branch that
    # is not identifiable from the data; PAML either objects or, on the Windows
    # build, aborts during setup without printing a reason. An unrooted tree is
    # written as a TRIFURCATION at the root, which is what the extra top-level
    # comma below produces. The biological rooting (Cricetidae outgroup) is
    # unchanged - it is simply not encoded as a root node.
    inner = f"({deo},{ger})" if ger else deo
    parts = [inner]
    if mur:
        parts.append(mur)
    if cri:
        parts.append(cri)
    if len(parts) == 1:
        return f"({parts[0]});"
    return "(" + ",".join(parts) + ");"


# ==============================================================================
# ALIGNMENT
# ==============================================================================

def load_cds(gene: str) -> dict[str, str]:
    gdir = cfg.ORTHO / gene
    if not gdir.exists():
        return {}
    out = {}
    for fa in gdir.glob("*.cds.fasta"):
        sp = fa.name.split(".")[0].replace("_", " ")
        try:
            s = str(next(SeqIO.parse(fa, "fasta")).seq).upper()
        except StopIteration:
            continue
        # MASK ambiguity rather than delete it. Deleting an N shifts every
        # downstream codon and manufactures internal stop codons, which the
        # frame QC would then blame on the annotation. (No sequence in the
        # current scaffold contains non-ACGT, so this is defensive.)
        s = "".join(c if c in "ACGT" else "N" for c in s)
        s = s[:len(s) - len(s) % 3]                    # trim to codons
        if len(s) >= MIN_CDS_LENGTH:
            out[sp] = s

    if USE_FIXED_PANEL:
        # All-or-nothing: a gene missing any panel member is dropped entirely,
        # rather than analysed on a smaller tree. A gene-specific topology is
        # exactly the confound this panel exists to remove.
        missing = [s for s in H2_PANEL if s not in out]
        if missing:
            log(f"  {gene}: dropped - panel incomplete, missing "
                f"{', '.join(missing)}", "WARN")
            return {}
        out = {s: out[s] for s in H2_PANEL}
    return out


# Reading-frame QC ------------------------------------------------------------
# Several CDSs in this scaffold came from Tier-2 BLAST HSP stitching, not from
# a native RefSeq CDS annotation. Stitching can join HSPs out of frame, which
# yields internal stop codons. codeml either aborts or - worse - runs and
# returns a confident dN/dS from a meaningless alignment. So frame is checked
# explicitly and failures are dropped with a named reason rather than silently
# passed through.

MAX_INTERNAL_STOPS = 0        # any internal stop disqualifies a sequence
MIN_CLEAN_FRACTION = 0.75     # drop the gene if >25% of species fail frame QC


def prepare_cds(sp: str, nt: str) -> tuple[str, str] | tuple[None, str]:
    """Return (trimmed_nt, peptide) or (None, reason)."""
    s = "".join(c if c in "ACGTacgt" else "N" for c in nt).upper()
    if len(s) < 60:
        return None, "too short (<20 codons)"
    trimmed = len(s) % 3
    if trimmed:
        s = s[:-trimmed]                       # trim a ragged 3' end
    pep = str(Seq(s).translate())
    if pep.endswith("*"):
        pep, s = pep[:-1], s[:-3]              # drop the terminal stop
    internal = pep.count("*")
    if internal > MAX_INTERNAL_STOPS:
        return None, f"{internal} internal stop codon(s) - frame is wrong"
    if pep.count("X") / max(len(pep), 1) > 0.10:
        return None, "greater than 10% ambiguous residues"
    return s, pep


def codon_align(gene: str, seqs: dict[str, str]) -> Path | None:
    """Translate -> align protein with MAFFT -> back-translate to codons.

    Codon-aware by construction. Aligning nucleotides directly would break the
    reading frame and make dN/dS meaningless - the single most common way to
    get a spurious selection result.
    """
    ALN.mkdir(parents=True, exist_ok=True)
    pep_in = ALN / f"{gene}.pep.fasta"
    pep_out = ALN / f"{gene}.pep.aln"
    cds_out = ALN / f"{gene}.codon.aln"
    if cds_out.exists():
        return cds_out

    clean: dict[str, str] = {}       # species -> frame-checked nt
    peps: dict[str, str] = {}
    for sp, s in seqs.items():
        nt, pep = prepare_cds(sp, s)
        if nt is None:
            log(f"  {gene}: DROP {sp} - {pep}", "WARN")
            continue
        clean[sp], peps[sp] = nt, pep

    if len(clean) < MIN_SPECIES_FOR_CODEML:
        log(f"  {gene}: only {len(clean)} species pass frame QC "
            f"(need {MIN_SPECIES_FOR_CODEML}) - skipping", "WARN")
        return None
    if len(clean) / max(len(seqs), 1) < MIN_CLEAN_FRACTION:
        log(f"  {gene}: {len(seqs) - len(clean)}/{len(seqs)} species fail "
            f"frame QC - scaffold quality too poor, skipping", "WARN")
        return None

    seqs = clean                      # back-translate against the trimmed nt
    with open(pep_in, "w") as fh:
        for sp, p in peps.items():
            fh.write(f">{sp.replace(' ', '_')}\n{p}\n")

    mafft = pc.tool("mafft")
    if mafft is None:
        log("  mafft not found", "WARN")
        return None
    try:
        with open(pep_out, "w") as fh:
            pc.run_tool(mafft, ["--auto", "--quiet", str(pep_in)],
                        stdout=fh, timeout=1800)
    except Exception as e:
        log(f"  mafft failed on {gene}: {e}", "WARN")
        return None

    # back-translate: walk the aligned protein, consuming codons
    aligned = {r.id: str(r.seq) for r in SeqIO.parse(pep_out, "fasta")}
    with open(cds_out, "w") as fh:
        for sp_key, apep in aligned.items():
            sp = sp_key.replace("_", " ")
            nt = seqs.get(sp, "")
            i, buf = 0, []
            for aa in apep:
                if aa == "-":
                    buf.append("---")
                else:
                    buf.append(nt[i:i + 3] if i + 3 <= len(nt) else "---")
                    i += 3
            fh.write(f">{sp_key}\n{''.join(buf)}\n")
    return cds_out


# The control file sets cleandata = 1, which makes codeml DELETE every column
# containing a gap or ambiguity before fitting. A ragged alignment can therefore
# leave codeml estimating dN/dS from a handful of codons - it will not complain,
# it will just return a confident, meaningless LRT. So count what codeml will
# actually keep, and refuse to run below a floor set in advance.
MIN_CLEAN_CODONS = 100


def clean_codon_count(aln: Path) -> int:
    """Codons surviving cleandata=1: no gap, no ambiguity, in any species."""
    recs = [str(r.seq) for r in SeqIO.parse(aln, "fasta")]
    if not recs:
        return 0
    L = min(len(s) for s in recs) // 3
    n = 0
    for i in range(L):
        cols = [s[i * 3:i * 3 + 3] for s in recs]
        if all(len(c) == 3 and set(c) <= set("ACGT") for c in cols):
            n += 1
    return n


def write_phylip(aln: Path, dest: Path):
    """PAML sequential PHYLIP. Validates before writing - PAML's complaints
    about malformed input are terse and easy to misattribute to the model."""
    recs = list(SeqIO.parse(aln, "fasta"))
    if not recs:
        return False
    lens = {len(r.seq) for r in recs}
    if len(lens) != 1:
        log(f"  {aln.stem}: ragged alignment, lengths {sorted(lens)}", "WARN")
        return False
    L = lens.pop()
    if L % 3:
        log(f"  {aln.stem}: length {L} is not a multiple of 3", "WARN")
        return False
    names = [r.id[:30] for r in recs]
    if len(set(names)) != len(names):
        log(f"  {aln.stem}: duplicate names after 30-char truncation", "WARN")
        return False
    with open(dest, "w") as fh:
        fh.write(f" {len(recs)} {L}\n")
        for r, nm in zip(recs, names):
            # PAML terminates the name at whitespace; two spaces is the
            # documented separator.
            fh.write(f"{nm}  {str(r.seq)}\n")
    return True


# ==============================================================================
# CODEML
# ==============================================================================

# A branch-site model A control file. Every option is stated explicitly rather
# than left to a default, because codeml's failure mode for an unsupported
# combination is to abort with exit status 1 and a message on STDOUT - which is
# invisible if stdout is captured and discarded.
#
# fix_alpha/alpha matter here: NSsites models estimate the omega distribution,
# so the gamma rate model must be switched OFF (alpha = 0, fix_alpha = 1).
# Leaving alpha free is an invalid combination with NSsites = 2.
CTL = """      seqfile = aln.phy
     treefile = tree.nwk
      outfile = out.txt
        noisy = 3
      verbose = 1
      runmode = 0
      seqtype = 1
    CodonFreq = 2
        clock = 0
       aaDist = 0
        model = 2
      NSsites = 2
        icode = 0
    fix_kappa = 0
        kappa = 2
    fix_omega = {fix_omega}
        omega = {omega}
    fix_alpha = 1
        alpha = 0
       Malpha = 0
        ncatG = 3
        getSE = 0
 RateAncestor = 0
   Small_Diff = .5e-6
    cleandata = 1
  fix_blength = 0
       method = 0
"""


def _codeml_diagnosis(out: str, err: str) -> str:
    """Turn codeml's message into something actionable."""
    blob = (out + "\n" + err).lower()
    for needle, msg in [
        ("perhaps in the tree", "tree file rejected - check the Newick and the #1 label"),
        ("error in tree", "tree file rejected by PAML's parser"),
        ("species not found", "a tip name in the tree is absent from the alignment"),
        ("sequence data file", "PAML could not read aln.phy"),
        ("sequence error", "PHYLIP parse error in aln.phy"),
        ("not a multiple of 3", "alignment length is not a multiple of 3"),
        ("stop codon", "a stop codon is present in the alignment"),
        ("ngene", "PHYLIP header line is malformed"),
    ]:
        if needle in blob:
            return msg
    return "no recognised PAML error string"


def run_codeml(gene: str, aln: Path, tree: str,
               verbose: bool = False) -> dict | None:
    codeml = pc.tool("codeml")
    if codeml is None:
        return None
    # MULTIPLE STARTING VALUES for the alternative model.
    #
    # Branch-site model A constrains omega2 >= 1, but the alternative was being
    # started at omega = 0.5 - BELOW the valid region. The optimiser then
    # clamps to the boundary and frequently never searches upward. Evidence
    # this was happening in the first run: 17/32 genes returned LRT exactly
    # 0.000, 18/32 had foreground omega sitting exactly at 1.0, and TIMP1
    # returned LRT = -0.001 - mathematically impossible for nested models, so
    # proof of non-convergence rather than absence of signal.
    #
    # PAML's documentation recommends multiple initial values because these
    # likelihood surfaces have local optima. The best (highest lnL) run is
    # kept, which is the maximum-likelihood estimate the test assumes it has.
    #
    # THIS IS NOT A FISHING EXPEDITION. It is applied identically to focal and
    # control genes, it can only move an under-converged alt model toward its
    # true optimum, and the pre-registered decision rule is COMPARATIVE -
    # focal must exceed all three control sets - so a uniform improvement in
    # convergence cannot favour the hypothesis.
    ALT_STARTS = (0.5, 1.5, 3.0)

    import shutil

    res = {}
    stages = [("null", 1, 1.0)] + [("alt", 0, w) for w in ALT_STARTS]
    for label, fix, om in stages:
        d = RUNS / gene / (label if label == "null" else f"alt_w{om}")
        # WIPE the run directory first. codeml leaves out.txt, rst, rst1, rub,
        # lnf and 2NG.* behind, and re-running into a dirty directory made the
        # results NON-DETERMINISTIC: two runs of identical code on identical
        # alignments returned different likelihoods (MS4A1 null moved +1.70,
        # PPIA alt moved +4.50, CD8A and CCR2 alt collapsed onto the null).
        # A likelihood fit must not depend on what a previous fit left on disk.
        if d.exists():
            shutil.rmtree(d, ignore_errors=True)
        d.mkdir(parents=True, exist_ok=True)
        if not write_phylip(aln, d / "aln.phy"):
            return None
        # PAML's tree file conventionally opens with "<ntaxa> <ntrees>".
        # Newer versions tolerate a bare Newick, older ones do not; the header
        # is harmless either way, so write it.
        ntax = len(list(SeqIO.parse(aln, "fasta")))
        deg = node_degrees(tree)
        if deg and (max(deg[:-1] or [0]) > 2 or deg[-1] > 3):
            log(f"  {gene}: tree has a polytomy (degrees {deg}) - codeml "
                "would abort with 'raise MAXNSONS'", "WARN")
            return None
        (d / "tree.nwk").write_text(f"{ntax} 1\n{tree}\n")
        (d / "codeml.ctl").write_text(CTL.format(fix_omega=fix, omega=om))
        try:
            pc.run_tool(codeml, ["codeml.ctl"], cwd=str(d),
                        capture_output=True, text=True, timeout=7200)
        except Exception as e:
            out = getattr(e, "stdout", "") or ""
            err = getattr(e, "stderr", "") or ""
            log(f"  codeml {label} failed for {gene}: "
                f"{_codeml_diagnosis(out, err)}", "WARN")
            if verbose or not getattr(run_codeml, "_shown", False):
                run_codeml._shown = True
                log(f"  --- codeml output ({gene}/{label}) ---", "WARN")
                for line in (out or err or str(e)).strip().splitlines()[-25:]:
                    log(f"    {line}", "WARN")
                log(f"  --- inputs are in {d} ---", "WARN")
            return None
        txt = (d / "out.txt").read_text() if (d / "out.txt").exists() else ""
        # Take the LAST lnL, not the first. PAML appends to outfile in some
        # configurations, so a stale block from an earlier invocation would
        # otherwise be read as the current result - silently returning the
        # previous run's answer.
        hits = re.findall(r"lnL\(.*?\):\s*(-?\d+\.\d+)", txt)
        if not hits:
            return None
        if len(hits) > 1:
            log(f"  {gene}/{label}: out.txt holds {len(hits)} lnL blocks "
                f"{hits} - taking the last. The directory was not clean.",
                "WARN")
        lnl = float(hits[-1])
        if label == "null":
            res["null"] = lnl
        else:
            # Keep the best (highest lnL) alternative run.
            if "alt" not in res or lnl > res["alt"]:
                res["alt"] = lnl
                w = re.search(r"foreground w\s+([\d.]+)\s+([\d.]+)\s+([\d.]+)",
                              txt)
                res["omega_fg"] = float(w.group(3)) if w else np.nan
                res["best_start"] = om
            res.setdefault("alt_lnL_by_start", {})[om] = lnl

    # A negative LRT means the alternative failed to reach the likelihood of
    # the null it contains - impossible at the true optimum, so a convergence
    # failure. Report it rather than clipping it silently to zero.
    if "alt" in res and "null" in res and res["alt"] < res["null"] - 1e-6:
        spread = res.get("alt_lnL_by_start", {})
        log(f"  {gene}: alt lnL ({res['alt']:.4f}) < null ({res['null']:.4f}) "
            f"even after {len(ALT_STARTS)} starting values {spread} - "
            "optimiser did not converge; treat this gene as uninformative "
            "rather than as evidence of no selection", "WARN")
    return res


# ==============================================================================
# DRIVER
# ==============================================================================

def gene_sets() -> dict[str, list[str]]:
    d = {"FOCAL": list(cfg.FOCAL_GENES)}
    # Secondary endpoint: E. Hassanat's a priori panel, three blocks tested
    # separately. Per her decision (20 Aug) WP2 is NOT re-run for the eight
    # genes absent from the ortholog scaffold - H1/H2 run on whatever is
    # available and the shortfall is reported, not repaired.
    for k, v in cfg.EXTENDED_PANEL.items():
        d[f"PANEL_{k}"] = list(v)
    d["PANEL_combined15"] = list(cfg.EXTENDED_15)
    for k, v in cfg.CONTROL_SETS.items():
        if v:
            d[k] = v
    return d


def panel_coverage_report():
    """What the extended panel can and cannot support at sequence level.

    Only 7 of the 15 panel genes exist in the ortholog scaffold, because the
    other 8 were never part of the original design and no orthologs were
    fetched. E.H. declined a WP2 re-run, so this shortfall is a reported
    limitation rather than something to fix.
    """
    log("\n" + "=" * 74)
    log("EXTENDED PANEL COVERAGE AT SEQUENCE LEVEL")
    log("=" * 74)
    testable = {}
    for blk, genes in cfg.EXTENDED_PANEL.items():
        rows = []
        for g in genes:
            n = len(load_cds(g))
            rows.append((g, n))
        ok = [g for g, n in rows if n >= MIN_SPECIES_FOR_CODEML]
        testable[blk] = ok
        log(f"\n  {blk}:")
        for g, n in rows:
            mark = "ok" if n >= MIN_SPECIES_FOR_CODEML else "NOT TESTABLE"
            log(f"    {g:10s} {n:2d} species  {mark}")
        log(f"    -> {len(ok)}/{len(genes)} testable "
            f"(module tests need >= 4)")
        if len(ok) < 4:
            log(f"    BLOCK NOT TESTABLE at sequence level - reported as a "
                f"data limitation, not a null result.", "WARN")
    return testable


def check() -> bool:
    log("=" * 74)
    log("WP4 H2 READINESS - branch-site selection")
    log("=" * 74)
    ok = True
    for t in ("mafft", "codeml"):
        p = pc.tool(t)
        log(f"  [{'OK     ' if p else 'MISSING'}] {t}  {p or ''}")
        ok &= p is not None
    if not ok:
        log("\n  Windows install (verified 19 Aug 2026):", "WARN")
        log("    codeml  https://github.com/abacus-gene/paml/releases/download/"
            "v4.10.10/paml-4.10.10-win-x86_64.tar.gz", "WARN")
        log("            The GitHub 'Source code' zip is source-ONLY - it has "
            "no bin\\ and no .exe. Use the win-x86_64 asset above.", "WARN")
        log("    mafft   https://mafft.cbrc.jp/alignment/software/"
            "mafft-7.526-win64-signed.zip", "WARN")
        log("            Ships mafft.bat, not mafft.exe - this is handled by "
            "platform_compat.run_tool().", "WARN")
        log("  Extract both under <drive>:\\Tools, or set ACOMYS_TOOL_DIRS.",
            "WARN")

    log("\nData:")
    n_ok = 0
    for name, genes in gene_sets().items():
        counts = [len(load_cds(g)) for g in genes]
        usable = sum(1 for c in counts if c >= MIN_SPECIES_FOR_CODEML)
        n_ok += usable
        log(f"  {name:20s} {usable}/{len(genes)} genes with "
            f">={MIN_SPECIES_FOR_CODEML} species")
    log(f"\n  {n_ok} genes total would enter the analysis")
    panel_coverage_report()
    if n_ok < 20:
        log("  Too few for a control-set comparison. Run "
            "02_build_ortholog_scaffold.py --all first.", "WARN")
        ok = False

    sp = sorted({s for g in cfg.FOCAL_GENES for s in load_cds(g)})
    if sp:
        log(f"\nTree (foreground = {FOREGROUND_CLADE}, marked #1):")
        log(f"  {build_tree(sp)}")
    log("\n" + "=" * 74)
    log(f"READY: {ok}")
    log("=" * 74)
    return ok


def align_all():
    log("=" * 74)
    log("CODON-AWARE ALIGNMENT (translate -> MAFFT -> back-translate)")
    log("=" * 74)
    thin = []
    for name, genes in gene_sets().items():
        for g in genes:
            seqs = load_cds(g)
            if len(seqs) < MIN_SPECIES_FOR_CODEML:
                continue
            p = codon_align(g, seqs)
            if not p:
                log(f"  {name:20s} {g:12s} {len(seqs)} species  FAILED")
                continue
            nc = clean_codon_count(p)
            flag = "" if nc >= MIN_CLEAN_CODONS else "  <-- TOO THIN"
            if flag:
                thin.append(g)
            log(f"  {name:20s} {g:12s} {len(seqs)} species  ok  "
                f"{nc:4d} clean codons{flag}")
    log(f"\nalignments in {ALN}")
    if thin:
        log(f"\n{len(thin)} alignment(s) below {MIN_CLEAN_CODONS} clean "
            f"codons and will be SKIPPED by --run:", "WARN")
        log(f"  {', '.join(thin)}", "WARN")
        log("  cleandata=1 deletes every gapped column, so these would give a "
            "confident LRT fitted on very little data.", "WARN")
    log("\nINSPECT A FEW BY EYE before running codeml. Automated codon")
    log("alignment fails silently on divergent chemokines, and a bad")
    log("alignment produces confident, wrong dN/dS.")


def run_all():
    log("=" * 74)
    log("CODEML branch-site model A - this is slow (minutes per gene)")
    log("=" * 74)
    rows = []
    # A gene can belong to several sets - all five PANEL_ECM_fibrosis genes are
    # also fibrosis_effector members. Without caching, codeml would run twice
    # per shared gene (wasted hours) and, worse, the same p-value would enter
    # the BH family twice, inflating the correction and double-counting the
    # gene's contribution. Run once per gene; report into every set it belongs
    # to.
    cache: dict[str, dict | None] = {}
    # Attrition must be reported per gene set. Every filter below is applied
    # identically to focal and control genes, but if it removes them at very
    # different rates the surviving comparison is no longer like-for-like.
    attrition: dict[str, dict[str, int]] = {}
    for name, genes in gene_sets().items():
        a = attrition.setdefault(
            name, {"total": len(genes), "no_aln": 0, "no_fg": 0,
                   "thin": 0, "codeml_failed": 0, "ran": 0})
        for g in genes:
            aln = ALN / f"{g}.codon.aln"
            if not aln.exists():
                a["no_aln"] += 1
                continue
            sp = [r.id.replace("_", " ") for r in SeqIO.parse(aln, "fasta")]
            if sum(1 for s in sp if cfg.CLADE.get(s) == FOREGROUND_CLADE) == 0:
                log(f"  {g}: no foreground species - skipped", "WARN")
                a["no_fg"] += 1
                continue
            nc = clean_codon_count(aln)
            if nc < MIN_CLEAN_CODONS:
                log(f"  {g}: only {nc} codons survive cleandata=1 "
                    f"(floor {MIN_CLEAN_CODONS}) - skipped", "WARN")
                a["thin"] += 1
                continue
            tree = build_tree(sp)
            if g in cache:
                r = cache[g]
            else:
                r = run_codeml(g, aln, tree)
                cache[g] = r
            if r is None:
                a["codeml_failed"] += 1
                continue
            a["ran"] += 1
            lrt = 2 * (r["alt"] - r["null"])
            p = stats.chi2.sf(max(lrt, 0), df=1)
            rows.append({"gene_set": name, "gene": g, "n_species": len(sp),
                         "lnL_null": r["null"], "lnL_alt": r["alt"],
                         "LRT": lrt, "omega_fg": r.get("omega_fg", np.nan),
                         "best_start": r.get("best_start", np.nan),
                         "converged": bool(r["alt"] >= r["null"] - 1e-6),
                         "p_raw": p})
            log(f"  {name:20s} {g:12s} LRT={lrt:7.3f}  p={p:.4g}")
    OUT.mkdir(parents=True, exist_ok=True)
    att = pd.DataFrame(attrition).T
    att["pct_ran"] = (100 * att["ran"] / att["total"]).round(1)
    att.to_csv(OUT / "attrition.csv")
    log("\n" + "=" * 74)
    log("GENE ATTRITION - read this before the results")
    log("=" * 74)
    log("\n" + att.to_string())
    rates = att.loc[att.total > 0, "pct_ran"]
    if len(rates) > 1 and (rates.max() - rates.min()) > 25:
        log("\nWARNING: gene sets survived filtering at very different rates "
            f"({rates.min():.0f}% to {rates.max():.0f}%).", "WARN")
        log("The Fisher comparison is then between sets of differing sequence "
            "quality, not differing selection. Report this alongside any "
            "positive result.", "WARN")

    if not rows:
        log("no results", "ERROR")
        return
    df = pd.DataFrame(rows)
    # BH across UNIQUE genes, then mapped back. Computing it on the row table
    # would count every shared gene once per set it belongs to and inflate n.
    uniq = df.drop_duplicates("gene")[["gene", "p_raw"]].copy()
    uniq["q_BH"] = _bh(uniq["p_raw"].values)
    df = df.merge(uniq[["gene", "q_BH"]], on="gene", how="left")
    log(f"\nBH computed over {len(uniq)} unique genes "
        f"({len(df)} set-memberships)")
    df.to_csv(OUT / "selection_results.csv", index=False)
    log(f"\nwrote {OUT / 'selection_results.csv'}")
    summarise()


def _bh(p):
    p = np.asarray(p, float); n = len(p); o = np.argsort(p)
    q = np.empty(n); prev = 1.0
    for rank, i in enumerate(o[::-1]):
        prev = min(prev, p[i] * n / (n - rank))
        q[i] = prev
    return q


def summarise():
    f = OUT / "selection_results.csv"
    if not f.exists():
        log("run --run first", "ERROR")
        return
    df = pd.read_csv(f)
    log("\n" + "=" * 74)
    log("H2 - BRANCH-SITE SELECTION ON THE ACOMYS LINEAGE")
    log("=" * 74)
    # omega_fg is only meaningful where the LRT is significant.
    #
    # When the MLE of omega2 sits at its lower bound of 1, the alternative
    # model collapses onto the null, the likelihood is flat in omega2, and the
    # reported value is an unidentifiable nuisance parameter - not an estimate.
    # Two byte-identical runs (same lnL to the last decimal) returned TUBB5
    # omega_fg = 4.75 and 58.56 respectively, both with LRT = 0.000 exactly.
    # Reporting a median omega over genes that mostly sit at the boundary
    # therefore averages numbers that mean nothing.
    df["omega_fg_reported"] = df["omega_fg"].where(df.q_BH < cfg.ALPHA)
    tab = (df.assign(sig=df.q_BH < cfg.ALPHA)
             .groupby("gene_set")
             .agg(n=("gene", "size"), n_sig=("sig", "sum"),
                  median_omega_sig=("omega_fg_reported", "median"),
                  n_at_boundary=("omega_fg", lambda s: int((s <= 1.0).sum()))))
    tab["pct_sig"] = (100 * tab.n_sig / tab.n).round(1)
    log("\n" + tab.to_string())
    log("\n  median_omega_sig is computed over SIGNIFICANT genes only; "
        "omega is unidentifiable where the LRT is 0.")
    log(f"  n_at_boundary counts genes with omega_fg <= 1 - the alternative "
        f"model has collapsed onto the null for those.")

    control_names = [k for k, v in cfg.CONTROL_SETS.items() if v]

    def fisher_block(name: str, controls: list[str]) -> tuple[int, int]:
        """Test one gene set against a named list of control sets."""
        if name not in tab.index:
            log(f"\n{name}: no genes entered the analysis - not testable "
                "at sequence level. Reported as a data limitation, NOT as a "
                "null result.", "WARN")
            return 0, 0
        fs, fn = int(tab.loc[name, "n_sig"]), int(tab.loc[name, "n"])
        usable = [c for c in controls if c in tab.index]
        if not usable:
            log(f"\n{name}: no usable control sets", "WARN")
            return 0, 0
        log(f"\nFisher's exact, {name} vs controls:")
        excluded = [c for c in control_names if c not in controls]
        if excluded:
            log(f"  (excluded for gene overlap: {', '.join(excluded)})")
        tests = []
        for c in usable:
            cs, cn = int(tab.loc[c, "n_sig"]), int(tab.loc[c, "n"])
            _, p = stats.fisher_exact([[fs, fn - fs], [cs, cn - cs]],
                                      alternative="greater")
            tests.append((c, p))
            log(f"  vs {c:22s} {fs}/{fn} vs {cs}/{cn}   p={p:.4g}")
        qs = _bh(np.array([p for _, p in tests]))
        for (c, _), q in zip(tests, qs):
            log(f"  q(BH) vs {c:20s} {q:.4g}")
        cleared = int((qs < cfg.ALPHA).sum())
        log(f"  -> cleared {cleared}/{len(tests)}")
        return cleared, len(tests)

    # ---- PRIMARY: the pre-registered focal module vs all three controls ----
    log("\n" + "=" * 74)
    log("PRIMARY ENDPOINT - pre-registered focal module")
    log("=" * 74)
    cleared, ntests = fisher_block("FOCAL", control_names)
    log("")
    log("PRE-SPECIFIED RULE: claim branch-specific selection ONLY IF the")
    log("focal module exceeds ALL THREE control sets at q<0.05.")
    log(f"RESULT: cleared {cleared}/{ntests} control sets -> "
        f"{'H2 SUPPORTED' if ntests and cleared == ntests else 'H2 NOT SUPPORTED'}")

    # ---- SECONDARY: E.H. panel, blocks tested separately ----
    log("\n" + "=" * 74)
    log("SECONDARY ENDPOINT - extended panel (E. Hassanat, 20 Aug 2026)")
    log("=" * 74)
    log("Blocks tested separately. Control sets sharing members with a block "
        "are excluded from that block's comparison (E.H. decision), because "
        "a Fisher test against a set containing half the block's own genes "
        "is self-referential.")
    for blk in cfg.EXTENDED_PANEL:
        fisher_block(f"PANEL_{blk}", cfg.controls_for(blk))
    log("\n  Combined 15-gene panel (supplementary):")
    fisher_block("PANEL_combined15",
                 [c for c in control_names
                  if c not in {"fibrosis_effector"}])
    log("=" * 74)


def debug_one(gene: str):
    """Run a single gene with everything shown. Use when --run fails."""
    aln = ALN / f"{gene}.codon.aln"
    if not aln.exists():
        log(f"no alignment for {gene} - run --align first", "ERROR")
        return
    recs = list(SeqIO.parse(aln, "fasta"))
    sp = [r.id.replace("_", " ") for r in recs]
    tree = build_tree(sp)

    log("=" * 74)
    log(f"DEBUG {gene}")
    log("=" * 74)
    log(f"\ncodeml : {pc.tool('codeml')}")
    log(f"species: {len(recs)}")
    for r in recs:
        log(f"  {r.id:26s} {len(r.seq):5d} nt")
    L = {len(r.seq) for r in recs}
    log(f"\nlengths: {sorted(L)}   multiple of 3: {all(x % 3 == 0 for x in L)}")
    log(f"clean codons (cleandata=1): {clean_codon_count(aln)}")
    log(f"\ntree:\n  {tree}")

    tips = set(re.findall(r"[A-Za-z]+_[A-Za-z]+", tree))
    ids = {r.id for r in recs}
    if tips != ids:
        log(f"\nMISMATCH between tree tips and alignment names:", "ERROR")
        log(f"  in tree only : {sorted(tips - ids)}", "ERROR")
        log(f"  in aln only  : {sorted(ids - tips)}", "ERROR")
    else:
        log("\ntree tips match alignment names exactly")
    log(f"'#1' marks: {tree.count('#1')}  (must be exactly 1)")

    d = node_degrees(tree)
    root, internal = d[-1], d[:-1]
    worst = max(internal) if internal else 0
    log(f"root degree: {root}  (3 = unrooted, which codeml wants)")
    log(f"max internal node degree: {worst}  (must be <= 2)")
    if worst > 2 or root > 3:
        log("POLYTOMY: codeml aborts with 'too many daughter nodes, raise "
            "MAXNSONS'. That message is printed after the pipe closes, so it "
            "does not appear in captured output.", "ERROR")

    log("\nrunning codeml with full output...\n")
    r = run_codeml(gene, aln, tree, verbose=True)
    if r:
        lrt = 2 * (r["alt"] - r["null"])
        log(f"\nOK  lnL null={r['null']:.4f}  alt={r['alt']:.4f}  "
            f"LRT={lrt:.4f}  p={stats.chi2.sf(max(lrt,0),df=1):.4g}")
    else:
        log(f"\nFAILED. Inspect {RUNS / gene / 'null'} - the .ctl, aln.phy "
            "and tree.nwk there are exactly what codeml was given.", "ERROR")


def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--check", action="store_true")
    ap.add_argument("--align", action="store_true")
    ap.add_argument("--run", action="store_true")
    ap.add_argument("--summarise", action="store_true")
    ap.add_argument("--debug", metavar="GENE",
                    help="run one gene with full codeml output")
    a = ap.parse_args()
    if a.check:
        sys.exit(0 if check() else 1)
    elif a.align:
        align_all()
    elif a.run:
        run_all()
    elif a.summarise:
        summarise()
    elif a.debug:
        debug_one(a.debug)
    else:
        ap.print_help()


if __name__ == "__main__":
    main()
