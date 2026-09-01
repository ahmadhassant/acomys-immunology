#!/usr/bin/env python3
"""
Trace resolve_gene() for one gene x species and print every decision point.

Written because three rounds of fixes to Tier 2 were made by inferring the
cause from a summary table. This shows the actual branch taken.

    python tests/debug_tier2.py                      # A. cahirinus / TGFB1
    python tests/debug_tier2.py --gene CCL2
    python tests/debug_tier2.py --species "Acomys dimidiatus"
"""

from __future__ import annotations

import argparse
import importlib.util
import os
import sys
from pathlib import Path
from shutil import which

HERE = Path(__file__).resolve().parent
CODE = HERE.parent
sys.path.insert(0, str(CODE))

import config as cfg

spec = importlib.util.spec_from_file_location(
    "wp2", CODE / "02_build_ortholog_scaffold.py")
wp2 = importlib.util.module_from_spec(spec)
sys.modules["wp2"] = wp2
spec.loader.exec_module(wp2)


def section(t):
    print(f"\n{'=' * 70}\n{t}\n{'=' * 70}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--gene", default="TGFB1")
    ap.add_argument("--species", default=cfg.FOCAL_SPECIES)
    args = ap.parse_args()

    sp = next((s for s in cfg.SPECIES_PANEL if s.name == args.species), None)
    if sp is None:
        print(f"unknown species: {args.species}")
        sys.exit(1)

    section(f"TRACE  {args.species}  /  {args.gene}")
    print(f"annotation flag : {sp.annotation}")
    print(f"assembly        : {sp.assembly}")
    print(f"usable          : {sp.usable}")
    print(f"is congener     : {sp.name == wp2.CONGENER_REFERENCE}")
    print(f"is focal        : {sp.name == cfg.FOCAL_SPECIES}")

    section("ENVIRONMENT")
    for t in ("makeblastdb", "blastn", "blastdbcmd"):
        print(f"  {t:12s} {which(t) or 'NOT FOUND'}")
    print(f"  ACOMYS_SCRATCH {os.environ.get('ACOMYS_SCRATCH') or 'NOT SET'}")
    print(f"  NCBI_API_KEY   {'set' if os.environ.get('NCBI_API_KEY') else 'not set'}")

    section("BRANCH 1 - native RefSeq")
    if sp.annotation == "refseq":
        d = wp2.fetch_native_refseq(args.gene, sp.taxid)
        print(f"  attempted. result: {'HIT' if d else 'None'}")
        if d:
            print(f"    accession {d['accession']}  cds={len(d['cds'])}")
    else:
        print(f"  SKIPPED (annotation='{sp.annotation}', not 'refseq')")

    section("BRANCH 2 - TSA (focal only)")
    tsa = cfg.REF / "transcriptome" / "acomys_15organ_TSA.fasta"
    print(f"  focal? {sp.name == cfg.FOCAL_SPECIES}   TSA exists? {tsa.exists()}")
    print("  SKIPPED" if not tsa.exists() else "  would attempt")

    section("BRANCH 3 - projected annotation fallback")
    if sp.annotation == "projected":
        d = wp2.fetch_native_refseq(args.gene, sp.taxid)
        print(f"  attempted (NCBI Gene lookup). result: {'HIT' if d else 'None'}")
        if d:
            print("  -> data found here, so TIER 2 WILL BE SKIPPED")
            print("     source='projected', confidence='low'")
    else:
        print(f"  SKIPPED (annotation='{sp.annotation}', not 'projected')")

    section("BRANCH 4 - Tier 2 homology")
    if sp.name == wp2.CONGENER_REFERENCE:
        print("  SKIPPED - this species IS the congener")
        return

    print(f"  target.assembly = {sp.assembly}")
    fna = wp2.unpack_genome(sp.assembly) if sp.assembly else None
    print(f"  unpack_genome -> {fna}")
    if fna is None:
        print("  STOP: no genomic FASTA")
        return

    db = wp2.genome_blast_db(sp.assembly)
    print(f"  genome_blast_db -> {db}")
    if db is None:
        print("  STOP: BLAST db could not be built")
        return
    print(f"  blast_db_exists -> {wp2.blast_db_exists(db)}")

    gdir = cfg.ORTHO / args.gene
    print(f"\n  query dir: {gdir}")
    print(f"  exists   : {gdir.exists()}")
    if gdir.exists():
        avail = sorted(p.name for p in gdir.glob("*.fasta"))
        print(f"  files    : {len(avail)}")
        for r in ("promoter", "utr5", "cds", "utr3"):
            q, qsp = wp2._pick_query(gdir, r, sp.name, wp2.CONGENER_REFERENCE)
            print(f"    {r:9s} query={qsp or 'NONE':24s} "
                  f"{q.name if q else '-'}")

    print("\n  running tier2_homology() ...")
    res = wp2.tier2_homology(args.gene, sp)
    if res is None:
        print("  RESULT: None  <-- this is why the species is excluded")
    else:
        print(f"  RESULT: {len(res['regions'])} region(s) recovered")
        for r, s in res["regions"].items():
            print(f"    {r:9s} {len(s):6d} bp   query={res['query_species'].get(r)}")


if __name__ == "__main__":
    main()
