#!/usr/bin/env python3
"""
================================================================================
TEST FIXTURE - synthetic ortholog scaffold with a KNOWN planted signal
================================================================================

Generates fake data in exactly the layout 02_build_ortholog_scaffold.py writes
and 03_kmer_phylo_controlled.py reads. Two purposes:

  1. Contract test. Proves WP2's output layout matches WP3's input expectation
     WITHOUT waiting for real downloads. If someone changes the filename
     convention in one script and not the other, this catches it immediately.

  2. Positive control. Plants a real divergence signal in the FOCAL module's
     3'UTR for Acomys and nowhere else. WP3 should detect it in the focal set,
     in the utr3 region, and NOT in the control sets or in CDS. If the pipeline
     cannot recover a signal we put there on purpose, it will not find a real
     one either.

Usage:
    python tests/make_synthetic_orthologs.py                 # with signal
    python tests/make_synthetic_orthologs.py --null          # no signal at all
    python tests/make_synthetic_orthologs.py --effect 0.25   # tune strength

    # then:
    cd code && python 03_kmer_phylo_controlled.py --arm residual --region utr3

EXPECTED OUTCOME (with signal, effect >= 0.15):
    utr3  : FOCAL residual_z clearly above all control sets
    cds   : no separation
    --null: no separation anywhere  <- guards against the pipeline inventing
                                       signal where none exists

WARNING: writes into data/reference/orthologs/. Use --outdir to keep it away
from real data, or delete the directory before a real run.
================================================================================
"""

from __future__ import annotations

import argparse
import random
import shutil
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
import config as cfg

log = cfg.log

BASES = "ACGT"
REGION_LENGTHS = {"promoter": 2000, "utr5": 300, "cds": 1200, "utr3": 900}

# Motifs planted in the focal module's 3'UTR only. AU-rich elements and a
# miRNA-like seed - the kind of thing the regulatory hypothesis predicts.
PLANTED_MOTIFS = ["ATTTATTTATTTA", "TATTTAT", "GCACTTA"]


def random_seq(n: int, rng: random.Random, gc: float = 0.45) -> str:
    w = [(1 - gc) / 2, gc / 2, gc / 2, (1 - gc) / 2]   # A C G T
    return "".join(rng.choices(BASES, weights=w, k=n))


def mutate(seq: str, rate: float, rng: random.Random) -> str:
    s = list(seq)
    for _ in range(int(len(s) * rate)):
        i = rng.randrange(len(s))
        s[i] = rng.choice(BASES)
    return "".join(s)


def plant(seq: str, motifs: list[str], n_copies: int, rng: random.Random) -> str:
    s = list(seq)
    for _ in range(n_copies):
        m = rng.choice(motifs)
        i = rng.randrange(0, max(1, len(s) - len(m)))
        s[i:i + len(m)] = list(m)
    return "".join(s)


def build(outdir: Path, effect: float, null: bool, seed: int):
    rng = random.Random(seed)
    species = cfg.USABLE_SPECIES

    gene_sets = {"FOCAL": list(cfg.FOCAL_GENES)}
    for name, genes in cfg.CONTROL_SETS.items():
        if genes:
            gene_sets[name] = genes

    n_written = 0
    for set_name, genes in gene_sets.items():
        is_focal = (set_name == "FOCAL") and not null
        for gene in genes:
            gdir = outdir / gene
            gdir.mkdir(parents=True, exist_ok=True)

            # ancestral sequence per region, shared across species
            ancestral = {r: random_seq(L, rng) for r, L in REGION_LENGTHS.items()}

            for sp in species:
                t = cfg.DIVERGENCE_MA.get(sp.name, 18.3)
                tag = sp.name.replace(" ", "_")

                for region, anc in ancestral.items():
                    # neutral divergence proportional to time; CDS constrained
                    rate = 0.0035 * t * (0.3 if region == "cds" else 1.0)
                    seq = mutate(anc, rate, rng)

                    # plant the signal: focal module, 3'UTR, Acomys only
                    if (is_focal and region == "utr3"
                            and sp.clade == "Deomyinae"):
                        n_copies = max(1, int(effect * 40))
                        seq = plant(seq, PLANTED_MOTIFS, n_copies, rng)

                    header = (f">{tag}|{gene}|{region}|SYNTHETIC"
                              f"|annotation=synthetic|len={len(seq)}")
                    (gdir / f"{tag}.{region}.fasta").write_text(
                        f"{header}\n{seq}\n")
                    n_written += 1

                # 'full' = utr5 + cds + utr3, matching the real script
                parts = []
                for region in ("utr5", "cds", "utr3"):
                    f = gdir / f"{tag}.{region}.fasta"
                    parts.append(f.read_text().split("\n")[1])
                full = "".join(parts)
                (gdir / f"{tag}.full.fasta").write_text(
                    f">{tag}|{gene}|full|SYNTHETIC|annotation=synthetic"
                    f"|len={len(full)}\n{full}\n")
                n_written += 1

    return n_written, sum(len(g) for g in gene_sets.values()), len(species)


def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--outdir", type=Path, default=cfg.ORTHO)
    ap.add_argument("--effect", type=float, default=0.30,
                    help="planted signal strength, 0-1 (default 0.30)")
    ap.add_argument("--null", action="store_true",
                    help="plant NO signal - pipeline should find nothing")
    ap.add_argument("--seed", type=int, default=cfg.RANDOM_SEED)
    ap.add_argument("--clean", action="store_true",
                    help="delete outdir first")
    args = ap.parse_args()

    if args.clean and args.outdir.exists():
        shutil.rmtree(args.outdir)
        log(f"cleaned {args.outdir}")

    args.outdir.mkdir(parents=True, exist_ok=True)

    log("=" * 74)
    log("SYNTHETIC ORTHOLOG FIXTURE")
    log(f"mode   : {'NULL (no signal)' if args.null else f'SIGNAL effect={args.effect}'}")
    log(f"outdir : {args.outdir}")
    log("=" * 74)

    n_files, n_genes, n_species = build(
        args.outdir, args.effect, args.null, args.seed)

    log(f"wrote {n_files} FASTA files "
        f"({n_genes} genes x {n_species} species x 5 regions)")
    log("")
    log("Expected result when you now run WP3:")
    if args.null:
        log("  NO separation between FOCAL and control sets, in any region.")
        log("  If WP3 reports q < 0.05 here, the test has FAILED and the")
        log("  pipeline is inventing signal.")
    else:
        log("  utr3 : FOCAL residual_z clearly above all three control sets")
        log("  cds  : no separation (signal was not planted there)")
    log("")
    log("  cd code && python 03_kmer_phylo_controlled.py --arm residual --region utr3")
    log("=" * 74)


if __name__ == "__main__":
    main()
