#!/usr/bin/env python3
"""
================================================================================
CHECK NCBI FOR NEW DATA RELEVANT TO THIS PROJECT
================================================================================

Queries NCBI directly rather than relying on literature search, which lags
submissions by weeks. Run it periodically; it reports what is new since a date.

    python code\\check_new_data.py                 # since 1 Aug 2026
    python code\\check_new_data.py --since 2026/01/01
    python code\\check_new_data.py --watch-only     # only the flagged gaps

WHAT IT WATCHES, AND WHY

  1. Lophuromys genome  -- HIGHEST VALUE if it appears.
     All three regenerative species in the panel are Acomys, so "regenerative"
     is confounded with "genus Acomys". Lophuromys is a deomyine OUTSIDE
     Acomys, confirmed regenerative (Riddell et al. 2025 PNAS). A nuclear
     assembly would break the confound outright. As of 18 Aug 2026: only
     mitochondrial genomes exist.

  2. Acomys spleen / immune tissue expression -- would revive the dropped arm.
     No public dataset has Acomys spleen expression. GSE168876 is kidney-only;
     PRJNA342864 pooled all 15 organs before sequencing.

  3. Acomys kidney / injury datasets -- more replicates would address the
     n = 3 per group limit that makes WP5 gene-level DE exploratory.

  4. Improved Acomys annotation -- a RefSeq annotation for GCA_029890205.1
     would remove the CAT-projection circularity, and full-length 3'UTRs
     would resolve the question H1 could not answer.

  5. Any new deomyine or gerbil assembly -- improves the phylogenetic control.

Requires: biopython. Uses NCBI_API_KEY from .env if present.
================================================================================
"""

from __future__ import annotations

import argparse
import sys
import time
from pathlib import Path

from Bio import Entrez

sys.path.insert(0, str(Path(__file__).resolve().parent))
import config as cfg

Entrez.email = cfg.ENTREZ_EMAIL
Entrez.api_key = cfg.ncbi_api_key()
DELAY = 0.12 if Entrez.api_key else 0.34
log = cfg.log

# taxids
ACOMYS_GENUS = 10068          # A. cahirinus; genus queries use the name
DEOMYINAE_GENERA = ["Acomys", "Lophuromys", "Deomys", "Uranomys"]
GERBIL_GENERA = ["Meriones", "Psammomys", "Gerbillus", "Rhombomys"]


def q(db: str, term: str, retmax: int = 40) -> list[str]:
    try:
        h = Entrez.esearch(db=db, term=term, retmax=retmax, sort="date")
        ids = Entrez.read(h)["IdList"]
        h.close()
        time.sleep(DELAY)
        return ids
    except Exception as e:
        log(f"    query failed ({db}): {e}", "WARN")
        return []


def summarise(db: str, ids: list[str], fields: tuple) -> list[dict]:
    if not ids:
        return []
    out = []
    try:
        h = Entrez.esummary(db=db, id=",".join(ids[:40]))
        recs = Entrez.read(h)
        h.close()
        time.sleep(DELAY)
        items = recs.get("DocumentSummarySet", {}).get("DocumentSummary", recs) \
            if isinstance(recs, dict) else recs
        for r in items:
            out.append({f: str(r.get(f, ""))[:110] for f in fields})
    except Exception as e:
        log(f"    summary failed ({db}): {e}", "WARN")
    return out


def section(title: str):
    log("")
    log("=" * 74)
    log(title)
    log("=" * 74)


def check_lophuromys():
    section("1. LOPHUROMYS GENOME  [highest value - breaks the Acomys confound]")
    ids = q("assembly", "Lophuromys[Organism]")
    if ids:
        log(f"  *** {len(ids)} ASSEMBLY RECORD(S) FOUND ***", "WARN")
        for r in summarise("assembly", ids,
                           ("AssemblyAccession", "AssemblyName",
                            "SpeciesName", "AssemblyStatus")):
            log(f"    {r['AssemblyAccession']}  {r['SpeciesName']}  "
                f"{r['AssemblyName']}  [{r['AssemblyStatus']}]")
        log("  -> If nuclear (not mitochondrial), ADD IT. This is the single")
        log("     most valuable species that could join the panel.")
    else:
        log("  none - only mitochondrial genomes exist (checked 18 Aug 2026)")

    for extra in ("Deomys", "Uranomys"):
        ids = q("assembly", f"{extra}[Organism]")
        log(f"  {extra}: {len(ids)} assembly record(s)")


def check_acomys_expression(since: str):
    section("2-3. ACOMYS EXPRESSION DATA  [spleen/immune, kidney, any new]")
    for label, term in [
        ("ALL Acomys GEO series",
         f'Acomys[Organism] AND ("{since}"[PDAT] : "3000"[PDAT])'),
        ("Acomys spleen / immune",
         'Acomys[Organism] AND (spleen[All Fields] OR immune[All Fields] '
         'OR macrophage[All Fields] OR monocyte[All Fields])'),
        ("Acomys kidney / renal",
         'Acomys[Organism] AND (kidney[All Fields] OR renal[All Fields] '
         'OR nephro[All Fields])'),
        ("Acomys single-cell",
         'Acomys[Organism] AND ("single cell"[All Fields] OR scRNA[All Fields])'),
    ]:
        ids = q("gds", term)
        log(f"\n  {label}: {len(ids)} hit(s)")
        for r in summarise("gds", ids, ("Accession", "title", "PDAT",
                                        "taxon", "n_samples")):
            log(f"    {r.get('Accession',''):12s} {r.get('PDAT',''):11s} "
                f"n={r.get('n_samples','?'):>3s}  {r.get('title','')[:70]}")


def check_annotation():
    section("4. ACOMYS ANNOTATION  [would resolve the H1 3'UTR limitation]")
    ids = q("assembly", "Acomys[Organism]")
    log(f"  {len(ids)} Acomys assembly record(s)")
    for r in summarise("assembly", ids,
                       ("AssemblyAccession", "SpeciesName", "AssemblyName",
                        "RefSeq_category", "AssemblyStatus")):
        flag = ""
        acc = r["AssemblyAccession"]
        if acc.startswith("GCF"):
            flag = "  <-- RefSeq-annotated"
        if acc.startswith("GCF") and "cahirinus" in r["SpeciesName"]:
            flag = "  <-- *** RefSeq for the FOCAL species: removes the " \
                   "CAT-projection circularity ***"
        log(f"    {acc:20s} {r['SpeciesName']:24s} {r['AssemblyName'][:22]:24s}"
            f"{flag}")


def check_comparators():
    section("5. NEW COMPARATOR GENOMES  [improve the phylogenetic control]")
    have = {s.assembly for s in cfg.SPECIES_PANEL if s.assembly}
    for genus in GERBIL_GENERA + ["Apodemus", "Mastomys", "Grammomys"]:
        ids = q("assembly", f"{genus}[Organism] AND latest[filter]")
        if not ids:
            continue
        rows = summarise("assembly", ids,
                         ("AssemblyAccession", "SpeciesName", "AssemblyStatus"))
        new = [r for r in rows if r["AssemblyAccession"] not in have]
        if new:
            log(f"\n  {genus}: {len(new)} not in the current panel")
            for r in new[:6]:
                log(f"    {r['AssemblyAccession']:20s} {r['SpeciesName']:26s}"
                    f" [{r['AssemblyStatus']}]")


def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--since", default="2026/08/01",
                    help="YYYY/MM/DD for date-filtered queries")
    ap.add_argument("--watch-only", action="store_true",
                    help="only the two flagged gaps (Lophuromys, spleen)")
    args = ap.parse_args()

    log("=" * 74)
    log("NCBI NEW-DATA CHECK")
    log(f"since        : {args.since}")
    log(f"API key      : {'set' if Entrez.api_key else 'NOT SET (slower)'}")
    log(f"current panel: {len(cfg.USABLE_SPECIES)} species")
    log("=" * 74)

    check_lophuromys()
    check_acomys_expression(args.since)
    if not args.watch_only:
        check_annotation()
        check_comparators()

    section("WHAT WOULD ACTUALLY CHANGE THE PROJECT")
    log("  Lophuromys nuclear genome  -> breaks 'regenerative = genus Acomys'")
    log("  Acomys spleen expression   -> revives the dropped spleen arm")
    log("  RefSeq annot. for A. cahirinus -> removes CAT circularity, gives")
    log("                                    full-length 3'UTRs, would let H1")
    log("                                    be re-tested properly")
    log("  More Acomys kidney replicates -> lifts WP5 above n=3 exploratory")
    log("")
    log("  Anything else is nice to have but does not change a conclusion.")
    log("=" * 74)


if __name__ == "__main__":
    main()
