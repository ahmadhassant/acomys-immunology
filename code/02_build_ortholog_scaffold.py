#!/usr/bin/env python3
"""
================================================================================
WP2 - ORTHOLOG SCAFFOLD CONSTRUCTION
Acomys cahirinus regeneration project
================================================================================

THE LINCHPIN. Everything in WP3, WP4 and WP5 depends on this producing a
defensible 1:1 ortholog set with correctly sliced sequence regions. Budget more
time than feels necessary.

WHAT IT PRODUCES
    data/reference/orthologs/<GENE>/<Species_name>.<region>.fasta
        region in {promoter, utr5, cds, utr3, full}
    data/reference/orthologs/ortholog_table.csv       one row per gene x species
    data/reference/orthologs/qc_report.csv            per-gene confidence
    data/reference/orthologs/random_background.txt    200 genes for the null

This is exactly the layout 03_kmer_phylo_controlled.py expects.

THE CIRCULARITY PROBLEM (read before running)
---------------------------------------------
A. cahirinus (GCA_029890205.1) is annotated by CAT projection FROM GRCm39.
Measuring divergence-from-mouse using mouse-projected gene models is circular:
projected UTR boundaries inherit mouse structure, which is the very thing WP3
measures.

This script therefore resolves each region under up to THREE independent
annotation sources and records which was used:

    (a) native      RefSeq annotation for the species                 BEST
    (b) congener    A. russatus GCF_903995435.1 boundaries, mapped     GOOD
                    (independently RefSeq-annotated, ~5 Ma divergence)
    (c) projected   CAT annotation from GRCm39                        SUSPECT
    (d) inferred    fixed window from CDS boundary                    FALLBACK

WP3 must be reported under (a)+(b) as primary and (c) as sensitivity. The
`annotation_source` column exists so that is auditable rather than assumed.

STRATEGY
    Tier 1  NCBI Gene symbol -> RefSeq mRNA + genomic coordinates. Fast, works
            for the 9 RefSeq-annotated species.
    Tier 2  Reciprocal best hits (BLAST/DIAMOND) against the assembly. Needed
            for A. cahirinus and anything Tier 1 misses.
    Tier 3  Search the de novo 15-organ TSA (PRJNA342864). Last resort for
            A. cahirinus; transcript evidence gives real UTRs, not projected.
    Then    Synteny check for watchlist genes; OrthoFinder cross-check.

Usage:
    python 02_build_ortholog_scaffold.py --check          # report readiness
    python 02_build_ortholog_scaffold.py --genes focal    # focal module only
    python 02_build_ortholog_scaffold.py --all
    python 02_build_ortholog_scaffold.py --all --annotation-source native

Requirements:
    pip install biopython pandas numpy tqdm
    conda install -c bioconda blast diamond orthofinder ncbi-datasets-cli
================================================================================
"""

from __future__ import annotations

import argparse

import hashlib
import json
import re
import subprocess
import sys
import time
import zipfile
from collections import defaultdict
from dataclasses import dataclass, asdict
from pathlib import Path
from shutil import which

import numpy as np
import pandas as pd
from Bio import Entrez, SeqIO
from Bio.Seq import Seq

# ---- shared config ------------------------------------------------------------
sys.path.insert(0, str(Path(__file__).resolve().parent))
import config as cfg
import platform_compat as pc

cfg.set_all_seeds()
Entrez.email = cfg.ENTREZ_EMAIL
import os
Entrez.api_key = cfg.ncbi_api_key()
DELAY = 0.12 if Entrez.api_key else 0.34

log = cfg.log

CONGENER_REFERENCE = "Acomys russatus"   # independently RefSeq-annotated

# ---- Tier 2 HSP filtering -----------------------------------------------------
# A gene locus occupies one genomic neighbourhood. HSPs scattered further than
# this are repeats or paralogs, not the same gene.
MAX_LOCUS_SPAN = 300_000      # bp around the best HSP
MAX_LENGTH_RATIO = 2.0        # recovered seq may not exceed 2x the query
MIN_PIDENT = 70.0             # below this, homology is not credible
# Adjacent exons routinely report a few bp of query overlap. Reject an HSP only
# when it substantially duplicates query already covered.
OVERLAP_TOLERANCE = 0.50      # fraction of the shorter HSP

# Biologically plausible region lengths. Anything outside these is reported,
# because a silently absurd UTR poisons WP3's primary endpoint.
PLAUSIBLE_LENGTH = {
    "promoter": (500, 5_000),
    "utr5":     (20, 10_000),
    "cds":      (150, 30_000),
    "utr3":     (50, 25_000),
}


# ==============================================================================
# RECORD TYPES
# ==============================================================================

@dataclass
class OrthologRecord:
    gene: str
    species: str
    clade: str
    tier: str                  # native_refseq | rbh | tsa | congener | failed
    annotation_source: str     # native | congener | projected | inferred | none
    accession: str = ""
    mrna_length: int = 0
    cds_length: int = 0
    utr5_length: int = 0
    utr3_length: int = 0
    promoter_length: int = 0
    has_complete_cds: bool = False
    n_paralogs_found: int = 0
    confidence: str = "unknown"   # high | medium | low | excluded
    flags: str = ""


# ==============================================================================
# TIER 1 - NATIVE REFSEQ
# ==============================================================================

def fetch_native_refseq(symbol: str, taxid: int) -> dict | None:
    """Resolve gene symbol -> best RefSeq mRNA, with CDS feature coordinates.

    GenBank format (not FASTA) so the CDS feature is available and UTRs can be
    sliced at real boundaries rather than guessed.
    """
    try:
        h = Entrez.esearch(
            db="gene",
            term=f'"{symbol}"[Gene Name] AND txid{taxid}[Organism] '
                 f'AND alive[prop]',
            retmax=10)
        ids = Entrez.read(h)["IdList"]
        h.close()
        time.sleep(DELAY)
        if not ids:
            return None

        n_hits = len(ids)
        h = Entrez.elink(dbfrom="gene", db="nuccore", id=ids[0],
                         linkname="gene_nuccore_refseqrna")
        links = Entrez.read(h)
        h.close()
        time.sleep(DELAY)
        if not links or not links[0].get("LinkSetDb"):
            return None

        nuc_ids = [l["Id"] for l in links[0]["LinkSetDb"][0]["Link"]]
        h = Entrez.efetch(db="nuccore", id=nuc_ids[0],
                          rettype="gb", retmode="text")
        rec = SeqIO.read(h, "genbank")
        h.close()
        time.sleep(DELAY)

        cds_feat = next((f for f in rec.features if f.type == "CDS"), None)
        if cds_feat is None:
            return None

        start = int(cds_feat.location.start)
        end = int(cds_feat.location.end)
        full = str(rec.seq).upper()

        return {
            "accession": rec.id,
            "full": full,
            "utr5": full[:start],
            "cds": full[start:end],
            "utr3": full[end:],
            "complete_cds": "partial" not in str(cds_feat.qualifiers).lower(),
            "n_paralogs": n_hits,
        }
    except Exception as e:
        _diagnose_entrez_error(e, symbol, taxid)
        return None


_ENTREZ_FATAL_REPORTED = False


def _diagnose_entrez_error(e: Exception, symbol: str, taxid: int):
    """Turn repeated per-gene Entrez failures into one actionable message.

    A bad NCBI_API_KEY makes EVERY request return HTTP 400. Logged per gene at
    DEBUG level that reads as hundreds of individual lookup failures, which
    looks like missing data rather than a broken credential - and quietly
    poisons the whole run.
    """
    global _ENTREZ_FATAL_REPORTED
    msg = str(e)

    if "400" in msg and not _ENTREZ_FATAL_REPORTED:
        _ENTREZ_FATAL_REPORTED = True
        key = os.environ.get("NCBI_API_KEY")
        log("", "ERROR")
        log("  " + "!" * 68, "ERROR")
        log("  NCBI returned HTTP 400 for an Entrez request.", "ERROR")
        if key:
            shown = key[:6] + "..." if len(key) > 8 else key
            log(f"  NCBI_API_KEY is set (starts '{shown}', len {len(key)}).",
                "ERROR")
            log("  An INVALID key makes every Entrez call return 400, so all", "ERROR")
            log("  Tier 1 lookups fail and every species looks like missing", "ERROR")
            log("  data. A real key is 36 hex characters.", "ERROR")
            log("", "ERROR")
            log("  FIX - either set a valid key or drop it entirely:", "ERROR")
            log("    unset NCBI_API_KEY && sed -i '/NCBI_API_KEY/d' ~/.bashrc",
                "ERROR")
            log("  The pipeline works without a key, just rate-limited to 3/s.",
                "ERROR")
        else:
            log("  No API key set, so this is likely rate limiting or a", "ERROR")
            log("  malformed query. Slow down or retry.", "ERROR")
        log("  " + "!" * 68, "ERROR")
        log("  STOPPING - continuing would record real genes as absent.",
            "ERROR")
        raise SystemExit(2)

    log(f"    native lookup failed ({symbol}, txid{taxid}): {msg}", "DEBUG")


def fetch_promoter(symbol: str, taxid: int, upstream: int) -> str | None:
    """Pull `upstream` bp 5' of the annotated TSS from the genomic record."""
    try:
        h = Entrez.esearch(db="gene",
                           term=f'"{symbol}"[Gene Name] AND txid{taxid}[Organism]',
                           retmax=1)
        ids = Entrez.read(h)["IdList"]
        h.close()
        time.sleep(DELAY)
        if not ids:
            return None

        h = Entrez.esummary(db="gene", id=ids[0])
        summ = Entrez.read(h)
        h.close()
        time.sleep(DELAY)

        info = summ["DocumentSummarySet"]["DocumentSummary"][0]
        loc = info.get("GenomicInfo")
        if not loc:
            return None
        gi = loc[0]
        chracc = gi["ChrAccVer"]
        start, stop = int(gi["ChrStart"]), int(gi["ChrStop"])

        if start <= stop:          # plus strand
            s, e, strand = max(0, start - upstream), start, 1
        else:                      # minus strand
            s, e, strand = start, start + upstream, 2

        h = Entrez.efetch(db="nuccore", id=chracc, rettype="fasta",
                          retmode="text", seq_start=s + 1, seq_stop=e,
                          strand=strand)
        rec = SeqIO.read(h, "fasta")
        h.close()
        time.sleep(DELAY)
        return str(rec.seq).upper()
    except Exception as e:
        _diagnose_entrez_error(e, symbol, taxid)
        return None


# ==============================================================================
# TIER 2 - RECIPROCAL BEST HITS
# ==============================================================================

def blast_safe(path: Path) -> str:
    """Space-free path for BLAST. Delegates to platform_compat.

    BLAST parses -in/-query/-db as space-separated file LISTS, so any path
    containing a space is read as several files. Windows uses an 8.3 short
    path; Linux stages a symlink. See platform_compat.safe_path.
    """
    return pc.safe_path(path, scratch=os.environ.get("ACOMYS_SCRATCH"))


def space_free_dir(preferred: Path) -> Path:
    """Writable directory with no spaces, for BLAST outputs."""
    return pc.space_free_dir(preferred)


def blast_db_exists(dbname: Path) -> bool:
    """True if a usable BLAST db is present.

    Genomes above ~1 GB produce MULTI-VOLUME databases: <db>.00.nin,
    <db>.01.nin ... plus a <db>.nal alias. Checking only for '<db>.nin' - as
    the first version of this did - reports 'missing' for a database that was
    built perfectly well, and silently rebuilds it every run.
    """
    p = str(dbname)
    return any(Path(p + s).exists() for s in (".nin", ".pin", ".nal", ".pal",
                                              ".00.nin", ".00.pin"))


def build_blast_db(fasta: Path, dbname: Path, dbtype: str = "nucl") -> bool:
    if blast_db_exists(dbname):
        return True
    if not pc.tool("makeblastdb"):
        log("  makeblastdb not found - install blast and symlink it", "WARN")
        return False

    size_gb = fasta.stat().st_size / 1e9
    log(f"    makeblastdb on {fasta.name} ({size_gb:.1f} GB) -> {dbname.name}")
    log(f"    ~1 min for a chromosome-scale assembly; cached after this")
    try:
        r = subprocess.run([pc.require("makeblastdb"), "-in", blast_safe(fasta),
                            "-dbtype", dbtype,
                            "-out", blast_safe(dbname), "-title", dbname.name,
                            "-max_file_sz", "3GB"],
                           check=True, capture_output=True, timeout=28800)
        if not blast_db_exists(dbname):
            log(f"    makeblastdb reported success but no db files appeared. "
                f"stdout: {r.stdout.decode()[:200]}", "ERROR")
            return False
        log("    db built")
        return True
    except subprocess.CalledProcessError as e:
        log(f"    makeblastdb FAILED: {e.stderr.decode()[:300]}", "ERROR")
    except subprocess.TimeoutExpired:
        log("    makeblastdb TIMED OUT after 8h. If the genome is on /mnt/c, "
            "move BLAST dbs to Linux-native disk: export ACOMYS_SCRATCH="
            "$HOME/acomys_scratch", "ERROR")
    return False


def unpack_genome(assembly: str) -> Path | None:
    """Unpack a `datasets` zip and return the genomic FASTA."""
    gdir = cfg.REF / "genomes"
    unpacked = gdir / assembly
    hits = list(unpacked.rglob("*_genomic.fna"))
    if hits:
        return hits[0]
    zipped = gdir / f"{assembly}.zip"
    if not zipped.exists():
        return None
    import zipfile
    log(f"    unpacking {assembly} ({zipped.stat().st_size:,} bytes)...")
    try:
        with zipfile.ZipFile(zipped) as zf:
            members = [n for n in zf.namelist() if n.endswith("_genomic.fna")]
            if not members:
                return None
            zf.extract(members[0], unpacked)
        hits = list(unpacked.rglob("*_genomic.fna"))
        return hits[0] if hits else None
    except Exception as e:
        log(f"    unpack failed: {e}", "WARN")
        return None


def genome_blast_db(assembly: str) -> Path | None:
    """Build (once) and return a nucleotide BLAST db for an assembly.

    Writes to ACOMYS_SCRATCH when set. On WSL2 this matters: makeblastdb on a
    2.3 GB genome writes several GB, and doing that across the /mnt/c 9p bridge
    is dramatically slower than on Linux-native disk.
    """
    fna = unpack_genome(assembly)
    if fna is None:
        log(f"    no genomic FASTA for {assembly}", "WARN")
        return None

    scratch = os.environ.get("ACOMYS_SCRATCH")
    dbdir = space_free_dir(Path(scratch) / "blastdb" if scratch else fna.parent)
    db = dbdir / f"{assembly}_db"

    if blast_db_exists(db):
        return db
    return db if build_blast_db(fna, db) else None


def tier2_homology(gene: str, target: cfg.Species,
                   congener: str = CONGENER_REFERENCE) -> dict | None:
    """Recover sequence for an UNANNOTATED genome by region-wise homology.

    WHY THIS EXISTS
    ---------------
    A. cahirinus (GCA_029890205.1) - the FOCAL species - has no NCBI Gene
    records. Its CAT annotation ships with the 2023 paper, not with Entrez, so
    the `datasets` package is genome-only and every symbol lookup fails. Same
    for A. dimidiatus (Ensembl-annotated). Without this function the focal
    species contributes nothing and WP3 cannot run at all.

    APPROACH
    --------
    Blast each ALREADY-EXTRACTED region of the congener (A. russatus, RefSeq,
    ~5 Ma diverged) separately against the target genome, and take the matching
    genomic segment as the target's version of that region:

        A_russatus.utr3.fasta  --blastn-->  A. cahirinus genome  -->  utr3

    Region-wise rather than whole-mRNA because the congener regions are already
    correctly segmented, so no CDS boundary has to be re-derived - and a wrong
    boundary here would silently corrupt WP3's primary endpoint.

    3'UTRs and promoters are typically single-exon and align contiguously. CDS
    is multi-exon; the top HSP is used, which is adequate for k-mer composition
    but is NOT a substitute for a real gene model.

    HONESTY
    -------
    Returns annotation_source='homology_congener', confidence 'medium' at best.
    02b_call_utrs.py supersedes this with read-based evidence once the RNA-seq
    is available.
    """
    if target.assembly is None:
        return None
    db = genome_blast_db(target.assembly)
    if db is None:
        return None

    cdir = cfg.ORTHO / gene
    out: dict = {"accession": f"homology:{target.assembly}", "regions": {},
                 "query_species": {}}

    for region in ("promoter", "utr5", "cds", "utr3"):
        q, qsp = _pick_query(cdir, region, target.name, congener)
        if q is None:
            continue
        seq = _blast_extract(q, db, region, verbose=True)
        if seq and len(seq) >= cfg.MIN_REGION_LENGTH:
            out["regions"][region] = seq
            out["query_species"][region] = qsp

    if not out["regions"]:
        return None
    out["full"] = "".join(out["regions"].get(r, "")
                          for r in ("utr5", "cds", "utr3"))
    return out


def _pick_query(gene_dir: Path, region: str, target_species: str,
                congener: str) -> tuple[Path | None, str | None]:
    """Choose the BLAST query for one region: congener first, else nearest.

    The congener is preferred, but RefSeq models are often partial - e.g.
    A. russatus has no annotated 3'UTR for TGFB1 while every other species
    does. Since 3'UTR is WP3's primary region, silently skipping it would drop
    the focal species from the primary analysis for that gene.

    So fall back to the phylogenetically NEAREST species that actually has the
    region. Which species supplied it is recorded per region, because a query
    from a more distant species yields a less reliable boundary and the reader
    is entitled to know which.
    """
    first = gene_dir / f"{congener.replace(' ', '_')}.{region}.fasta"
    if first.exists():
        return first, congener

    t_div = cfg.DIVERGENCE_MA.get(target_species, 0.0)
    candidates = []
    for fa in gene_dir.glob(f"*.{region}.fasta"):
        sp = fa.name.split(".")[0].replace("_", " ")
        if sp == target_species:
            continue
        d = cfg.DIVERGENCE_MA.get(sp)
        if d is None:
            continue
        candidates.append((abs(d - t_div), sp, fa))

    if not candidates:
        return None, None
    candidates.sort(key=lambda c: c[0])
    _, sp, fa = candidates[0]
    return fa, sp


def _blast_extract(query: Path, db: Path, region: str,
                   evalue: float = 1e-5, verbose: bool = False) -> str | None:
    """blastn one region and return the aligned SUBJECT sequence.

    Takes `sseq` straight from the BLAST tabular output rather than fetching
    coordinates and then calling blastdbcmd. The blastdbcmd route needs the
    database to have been created with -parse_seqids; without it, `-entry`
    silently returns nothing, which is how this function previously produced
    no sequence at all despite BLAST finding good hits.

    `sseq` is the aligned subject with gap characters, which are stripped. For
    k-mer composition that is exactly what we want: the homologous segment,
    without introducing flanking sequence that has no counterpart in the query.

    Multiple HSPs on the same subject and strand are concatenated in query
    order, which recovers multi-exon CDS reasonably. It is NOT a gene model -
    confidence stays capped at 'medium'.
    """
    fields = ("6 sseqid sstart send pident length bitscore sstrand "
              "qstart qend sseq")
    try:
        res = subprocess.run(
            [pc.require("blastn"), "-query", blast_safe(query), "-db", blast_safe(db),
             "-evalue", str(evalue), "-max_target_seqs", "5",
             "-outfmt", fields],
            check=True, capture_output=True, timeout=3600).stdout.decode()
    except (subprocess.CalledProcessError, subprocess.TimeoutExpired,
            FileNotFoundError) as e:
        if verbose:
            log(f"      blastn failed for {region}: {e}", "WARN")
        return None

    rows = [l.split("\t") for l in res.strip().splitlines() if l]
    if not rows:
        if verbose:
            log(f"      no BLAST hits for {region}", "WARN")
        return None

    rows = [r for r in rows if len(r) >= 10 and float(r[3]) >= MIN_PIDENT]
    if not rows:
        if verbose:
            log(f"      {region}: hits all below {MIN_PIDENT}% identity", "WARN")
        return None

    # ---- anchor on the single best HSP ----
    rows.sort(key=lambda r: float(r[5]), reverse=True)
    anchor = rows[0]
    a_seqid, a_strand = anchor[0], anchor[6]
    a_mid = (int(anchor[1]) + int(anchor[2])) / 2

    # ---- keep only HSPs at the SAME LOCUS ----
    # Without this, a repeat-rich query (promoters and UTRs especially) hits
    # hundreds of positions across a chromosome and concatenating them yields
    # nonsense: a 2 kb promoter query produced 463 kb of "sequence" before this
    # filter existed. Restrict to one genomic neighbourhood around the anchor.
    local = []
    for r in rows:
        if r[0] != a_seqid or r[6] != a_strand:
            continue
        mid = (int(r[1]) + int(r[2])) / 2
        if abs(mid - a_mid) <= MAX_LOCUS_SPAN:
            local.append(r)

    # ---- enforce colinearity AND drop query-overlapping HSPs ----
    # Overlap matters as much as colinearity: two HSPs covering the same part
    # of the query are alternative alignments of one segment, not two exons.
    # Concatenating both duplicates sequence - CCR2 came out at 1.81x its true
    # CDS length that way, which is small enough to look plausible.
    # Use BLAST's own qstart/qend rather than inferring the query span from the
    # subject length, and tolerate the few bp of overlap that adjacent exons
    # normally report. Inferring the span and rejecting ANY overlap discarded
    # genuine exons - TGFB1/TGFB3 CDS fell to 0.65x of the query.
    local.sort(key=lambda r: (int(r[7]), -float(r[5])))   # qstart, best first
    colinear, last_s, covered = [], None, []
    for r in local:
        q_start, q_end = int(r[7]), int(r[8])
        if q_start > q_end:
            q_start, q_end = q_end, q_start
        qlen_hsp = q_end - q_start + 1

        redundant = False
        for c0, c1 in covered:
            ov = min(q_end, c1) - max(q_start, c0) + 1
            if ov > 0 and ov > OVERLAP_TOLERANCE * min(qlen_hsp, c1 - c0 + 1):
                redundant = True
                break
        if redundant:
            continue

        s_start = int(r[1])
        if last_s is not None:
            forward = a_strand != "minus"
            if (forward and s_start < last_s) or (not forward and s_start > last_s):
                continue                                 # breaks colinearity
        colinear.append(r)
        covered.append((q_start, q_end))
        last_s = s_start
    if not colinear:
        colinear = [anchor]

    seq = "".join(h[9].replace("-", "") for h in colinear).upper()
    seq = "".join(c for c in seq if c in "ACGTN")

    # ---- length sanity against the query ----
    qlen = _fasta_seq_len(query)
    if qlen and len(seq) > MAX_LENGTH_RATIO * qlen:
        # Still implausible - fall back to the anchor HSP alone.
        seq = "".join(c for c in anchor[9].replace("-", "").upper()
                      if c in "ACGTN")
        if verbose:
            log(f"      {region}: stitched length exceeded "
                f"{MAX_LENGTH_RATIO}x query ({qlen} bp) - using anchor HSP "
                f"only ({len(seq)} bp)", "WARN")

    if verbose:
        log(f"      {region}: {len(colinear)}/{len(rows)} HSP(s) kept on "
            f"{a_seqid} ({a_strand}), {len(seq)} bp "
            f"[query {qlen} bp, {float(anchor[3]):.1f}% id]")
    return seq or None


def _fasta_seq_len(path: Path) -> int:
    try:
        return sum(len(l.strip()) for l in path.read_text().splitlines()
                   if l and not l.startswith(">"))
    except Exception:
        return 0


def reciprocal_best_hit(query_fasta: Path, target_db: Path,
                        target_fasta: Path, query_db: Path,
                        evalue: float = 1e-20) -> dict[str, str]:
    """Forward and reverse BLAST; keep only mutually-best pairs.

    Reciprocity is what separates an ortholog from 'the most similar paralog'.
    Single-direction best hits are the classic way to build a wrong ortholog set.
    """
    if not pc.tool("blastn"):
        log("  blastn not found", "WARN")
        return {}

    def run(q, db):
        try:
            out = subprocess.run(
                [pc.require("blastn"), "-query", blast_safe(q), "-db", blast_safe(db),
                 "-evalue", str(evalue), "-max_target_seqs", "5",
                 "-outfmt", "6 qseqid sseqid pident length evalue bitscore"],
                check=True, capture_output=True, timeout=3600).stdout.decode()
        except subprocess.CalledProcessError:
            return {}
        best, seen = {}, {}
        for line in out.strip().splitlines():
            if not line:
                continue
            p = line.split("\t")
            q_id, s_id, bits = p[0], p[1], float(p[5])
            if q_id not in seen or bits > seen[q_id]:
                seen[q_id] = bits
                best[q_id] = s_id
        return best

    fwd = run(query_fasta, target_db)
    rev = run(target_fasta, query_db)
    return {q: s for q, s in fwd.items() if rev.get(s) == q}


# ==============================================================================
# TIER 3 - DE NOVO TSA SEARCH (Acomys only)
# ==============================================================================

def search_tsa(symbol: str, tsa_fasta: Path, query_seq: str,
               evalue: float = 1e-30) -> dict | None:
    """Find the gene in the 15-organ de novo assembly.

    Worth the effort for A. cahirinus: TSA contigs are transcript-derived, so
    the UTRs are real observed sequence rather than projected from mouse. This
    is the cleanest route around the circularity problem.
    """
    if not tsa_fasta.exists():
        return None
    db = tsa_fasta.parent / "acomys_tsa_db"
    if not build_blast_db(tsa_fasta, db):
        return None

    tmp_q = tsa_fasta.parent / f"_q_{symbol}.fasta"
    tmp_q.write_text(f">{symbol}\n{query_seq}\n")
    try:
        out = subprocess.run(
            [pc.require("blastn"), "-query", blast_safe(tmp_q), "-db", blast_safe(db),
             "-evalue", str(evalue), "-max_target_seqs", "3",
             "-outfmt", "6 sseqid pident length evalue bitscore"],
            check=True, capture_output=True, timeout=1800).stdout.decode()
    except subprocess.CalledProcessError:
        return None
    finally:
        tmp_q.unlink(missing_ok=True)

    rows = [l.split("\t") for l in out.strip().splitlines() if l]
    if not rows:
        return None
    hit_id = rows[0][0]
    for rec in SeqIO.parse(tsa_fasta, "fasta"):
        if rec.id == hit_id or rec.id.split()[0] == hit_id:
            return {"accession": hit_id, "full": str(rec.seq).upper(),
                    "pident": float(rows[0][1]), "n_hits": len(rows)}
    return None


def predict_orf(seq: str) -> tuple[int, int] | None:
    """Longest ORF, all three frames, both strands. Used to derive UTR
    boundaries from a TSA contig that carries no CDS annotation."""
    best = None
    for strand, s in ((1, seq), (-1, str(Seq(seq).reverse_complement()))):
        for frame in range(3):
            trimmed = s[frame:len(s) - ((len(s) - frame) % 3)]
            if len(trimmed) < 3:
                continue
            prot = str(Seq(trimmed).translate())
            for m in re.finditer(r"M[^*]*\*", prot):
                length = m.end() - m.start()
                if best is None or length > best[0]:
                    nt_s = frame + m.start() * 3
                    nt_e = frame + m.end() * 3
                    best = (length, strand, nt_s, nt_e)
    if best is None or best[0] < 50:
        return None
    _, strand, s0, e0 = best
    if strand == -1:
        n = len(seq)
        return (n - e0, n - s0)
    return (s0, e0)


# ==============================================================================
# REGION SLICING
# ==============================================================================

def slice_regions(data: dict, promoter: str | None,
                  annotation_source: str) -> dict[str, str]:
    """Assemble the five region sequences, filling gaps by fixed windows."""
    full = data.get("full", "")
    cds = data.get("cds", "")
    utr5 = data.get("utr5", "")
    utr3 = data.get("utr3", "")

    # No CDS annotation (TSA route) -> predict the ORF.
    if not cds and full:
        orf = predict_orf(full)
        if orf:
            s, e = orf
            cds, utr5, utr3 = full[s:e], full[:s], full[e:]

    # Short/absent UTRs -> fixed window from the CDS boundary, flagged inferred.
    L = cfg.FALLBACK_UTR_LENGTH
    if len(utr5) < cfg.MIN_REGION_LENGTH and full and cds:
        i = full.find(cds)
        if i > 0:
            utr5 = full[max(0, i - L):i]
    if len(utr3) < cfg.MIN_REGION_LENGTH and full and cds:
        i = full.find(cds)
        if i >= 0:
            j = i + len(cds)
            utr3 = full[j:j + L]

    return {"promoter": promoter or "", "utr5": utr5, "cds": cds,
            "utr3": utr3, "full": full}


def write_regions(gene: str, species: str, regions: dict[str, str],
                  accession: str, annotation_source: str) -> dict[str, int]:
    gdir = cfg.ORTHO / gene
    gdir.mkdir(parents=True, exist_ok=True)
    tag = species.replace(" ", "_")
    lengths = {}
    for region, seq in regions.items():
        lengths[region] = len(seq)
        if len(seq) < cfg.MIN_REGION_LENGTH:
            continue
        header = (f">{tag}|{gene}|{region}|{accession}"
                  f"|annotation={annotation_source}|len={len(seq)}")
        (gdir / f"{tag}.{region}.fasta").write_text(f"{header}\n{seq}\n")
    return lengths


# ==============================================================================
# PER-GENE RESOLUTION
# ==============================================================================

def resolve_gene(gene: str, species: cfg.Species,
                 tsa_fasta: Path | None,
                 congener_cache: dict) -> OrthologRecord:
    rec = OrthologRecord(gene=gene, species=species.name, clade=species.clade,
                         tier="failed", annotation_source="none")
    flags = []

    # ---- Tier 1: native RefSeq ----
    data, source, tier = None, "none", "failed"
    if species.annotation == "refseq":
        for cand in cfg.symbol_candidates(gene):
            data = fetch_native_refseq(cand, species.taxid)
            if data:
                source, tier = "native", "native_refseq"
                if cand != gene:
                    flags.append(f"resolved_via_alias={cand}")
                break

    # ---- Tier 3 before Tier 1-projected, for Acomys cahirinus ----
    # Deliberate ordering: transcript evidence beats mouse-projected models.
    if data is None and species.name == cfg.FOCAL_SPECIES and tsa_fasta:
        query = congener_cache.get(gene) or ""
        if query:
            hit = search_tsa(gene, tsa_fasta, query)
            if hit:
                data = hit
                source, tier = "congener", "tsa"
                flags.append(f"tsa_pident={hit.get('pident', 0):.1f}")
                if hit.get("n_hits", 1) > 1:
                    flags.append(f"tsa_multihit={hit['n_hits']}")

    # ---- Tier 1 fallback: projected annotation (suspect) ----
    if data is None and species.annotation == "projected":
        data = fetch_native_refseq(gene, species.taxid)
        if data:
            source, tier = "projected", "native_refseq"
            flags.append("CIRCULARITY_RISK_mouse_projected")

    # ---- Tier 2: homology from the congener against an unannotated genome ----
    # Required for A. cahirinus (focal) and A. dimidiatus, whose annotations
    # are not in NCBI Gene. Without this the focal species contributes nothing.
    homology = None
    if data is None and species.name != CONGENER_REFERENCE:
        homology = tier2_homology(gene, species)
        if homology:
            source, tier = "homology_congener", "rbh"
            flags.append(f"homology_from_{CONGENER_REFERENCE.replace(' ', '_')}")
            data = {"accession": homology["accession"],
                    "full": homology.get("full", ""),
                    "cds": homology["regions"].get("cds", ""),
                    "utr5": homology["regions"].get("utr5", ""),
                    "utr3": homology["regions"].get("utr3", ""),
                    "complete_cds": False, "n_paralogs": 1}

    if data is None:
        rec.confidence = "excluded"
        rec.flags = ";".join(flags + ["no_sequence_found"])
        return rec

    if homology:
        promoter = homology["regions"].get("promoter")
    else:
        promoter = None
        for cand in cfg.symbol_candidates(gene):
            promoter = fetch_promoter(cand, species.taxid,
                                      cfg.PROMOTER_UPSTREAM_BP)
            if promoter:
                break
    if not promoter:
        flags.append("no_promoter")

    regions = slice_regions(data, promoter, source)
    lengths = write_regions(gene, species.name, regions,
                            data.get("accession", "NA"), source)

    if lengths.get("utr3", 0) and not data.get("utr3"):
        flags.append("utr3_inferred_window")
    if lengths.get("utr5", 0) and not data.get("utr5"):
        flags.append("utr5_inferred_window")

    # ---- flag biologically implausible region lengths ----
    for region, (lo, hi) in PLAUSIBLE_LENGTH.items():
        L = lengths.get(region, 0)
        if L and not (lo <= L <= hi):
            flags.append(f"IMPLAUSIBLE_{region}_len={L}")

    n_par = data.get("n_paralogs", 1)
    if gene in cfg.ORTHOLOGY_WATCHLIST:
        flags.append("WATCHLIST_verify_synteny")
    if n_par > 3:
        flags.append(f"many_gene_hits={n_par}")

    if any(f.startswith("IMPLAUSIBLE_") for f in flags):
        conf = "low"
    elif source == "native" and lengths.get("cds", 0) > 150 and n_par <= 2:
        conf = "high"
    elif source in ("native", "congener", "homology_congener"):
        conf = "medium"
    else:
        conf = "low"

    rec.tier = tier
    rec.annotation_source = source
    rec.accession = data.get("accession", "NA")
    rec.mrna_length = lengths.get("full", 0)
    rec.cds_length = lengths.get("cds", 0)
    rec.utr5_length = lengths.get("utr5", 0)
    rec.utr3_length = lengths.get("utr3", 0)
    rec.promoter_length = lengths.get("promoter", 0)
    rec.has_complete_cds = bool(data.get("complete_cds", False))
    rec.n_paralogs_found = n_par
    rec.confidence = conf
    rec.flags = ";".join(flags)
    return rec


# ==============================================================================
# RANDOM BACKGROUND
# ==============================================================================

def draw_random_background(n: int = cfg.N_RANDOM_BACKGROUND) -> list[str]:
    """Draw n mouse gene symbols at random as the empirical null gene set.

    Without this the residual test has no background distribution and 'unusual'
    has no referent.
    """
    out_file = cfg.ORTHO / "random_background.txt"
    if out_file.exists():
        genes = [l.strip() for l in out_file.read_text().splitlines() if l.strip()]
        log(f"  cached random background: {len(genes)} genes")
        return genes

    log(f"  drawing {n} random 1:1 background genes ...")
    try:
        h = Entrez.esearch(db="gene",
                           term="txid10090[Organism] AND alive[prop] "
                                "AND genetype protein coding[prop]",
                           retmax=6000, retstart=0)
        ids = Entrez.read(h)["IdList"]
        h.close()
        time.sleep(DELAY)
        rng = np.random.default_rng(cfg.RANDOM_SEED)
        pick = rng.choice(ids, size=min(n * 2, len(ids)), replace=False)

        symbols = []
        for i in range(0, len(pick), 200):
            batch = list(pick[i:i + 200])
            h = Entrez.esummary(db="gene", id=",".join(batch))
            summ = Entrez.read(h)
            h.close()
            time.sleep(DELAY)
            for d in summ["DocumentSummarySet"]["DocumentSummary"]:
                s = d.get("Name", "").upper()
                if s and s not in cfg.all_genes():
                    symbols.append(s)
            if len(symbols) >= n:
                break

        symbols = sorted(set(symbols))[:n]
        cfg.ORTHO.mkdir(parents=True, exist_ok=True)
        out_file.write_text("\n".join(symbols))
        log(f"  drew {len(symbols)} background genes")
        return symbols
    except Exception as e:
        log(f"  background draw failed: {e}", "WARN")
        return []


# ==============================================================================
# READINESS CHECK
# ==============================================================================

def check_readiness() -> bool:
    log("=" * 74)
    log("WP2 READINESS CHECK")
    log("=" * 74)
    ok = True

    log("\nTools:")
    win_hint = {
        "makeblastdb": "NCBI BLAST+ win64 installer",
        "blastn":      "NCBI BLAST+ win64 installer",
        "datasets":    "conda install -c conda-forge ncbi-datasets-cli",
        "diamond":     "not available natively on Windows (optional)",
        "orthofinder": "not available natively on Windows (optional)",
    }
    for tool, pkg in [("makeblastdb", "blast"), ("blastn", "blast"),
                      ("datasets", "ncbi-datasets-cli"),
                      ("diamond", "diamond"), ("orthofinder", "orthofinder")]:
        present = pc.tool(tool) is not None
        need = tool in ("makeblastdb", "blastn")
        mark = "OK " if present else ("MISSING" if need else "optional")
        hint = win_hint[tool] if pc.IS_WINDOWS else f"conda install -c bioconda {pkg}"
        log(f"  [{mark:8s}] {tool:14s} ({hint})")
        if need and not present:
            ok = False

    log("\nData:")
    tsa = cfg.REF / "transcriptome" / "acomys_15organ_TSA.fasta"
    genomes = list((cfg.REF / "genomes").glob("*.zip")) if (cfg.REF / "genomes").exists() else []
    checks = [
        ("Acomys transcript evidence (optional - congener annotation is used "
         "instead)", tsa.exists(),
         "not required; A. russatus/A. dimidiatus supply UTR boundaries", False),
        (f"Genome archives ({len(genomes)} found)", len(genomes) > 0,
         "python 01_fetch_data.py --step genomes", False),
        ("NCBI_API_KEY set", Entrez.api_key is not None,
         "export NCBI_API_KEY=... (3x faster)", False),
    ]
    for name, present, fix, required in checks:
        mark = "OK " if present else ("MISSING" if required else "optional")
        log(f"  [{mark:8s}] {name}")
        if not present:
            log(f"              -> {fix}")
        if required and not present:
            ok = False

    log("\nPanel:")
    n = len(cfg.USABLE_SPECIES)
    log(f"  {n} usable species (minimum {cfg.MIN_SPECIES_REQUIRED})")
    log(f"  GATE: {'PASS' if n >= cfg.MIN_SPECIES_REQUIRED else 'FAIL'}")
    if n < cfg.MIN_SPECIES_REQUIRED:
        ok = False

    log("\nCarried design gaps:")
    for key in cfg.DESIGN_GAPS:
        log(f"  * {key}")

    log("\n" + "=" * 74)
    log(f"READY: {ok}")
    log("=" * 74)
    return ok


# ==============================================================================
# MAIN
# ==============================================================================

def tier2_preflight(species_list: list) -> dict:
    """Check ONCE whether homology recovery can work, and say so loudly.

    Species with no NCBI annotation depend entirely on Tier 2. If the BLAST
    tooling is missing, they are silently excluded and the focal species can
    vanish from the analysis without any obvious error - which is exactly what
    happened on the first three runs of this script.
    """
    needs = [s for s in species_list if s.annotation not in ("refseq",)]
    status = {"needed_by": [s.name for s in needs], "ok": True, "reasons": []}
    if not needs:
        return status

    log("")
    log("  TIER 2 PREFLIGHT - homology recovery for unannotated genomes")
    log(f"    required by: {', '.join(s.name for s in needs)}")

    for tool in ("makeblastdb", "blastn", "blastdbcmd"):
        found = pc.tool(tool)   # searches Program Files on Windows, not just PATH
        log(f"    [{'OK     ' if found else 'MISSING'}] {tool}"
            f"{'  ' + found if found else ''}")
        if not found:
            status["ok"] = False
            status["reasons"].append(f"{tool} not on PATH")

    scratch = os.environ.get("ACOMYS_SCRATCH")
    log(f"    ACOMYS_SCRATCH: {scratch or 'NOT SET (dbs go on /mnt/c - slow)'}")

    for s in needs:
        fna = unpack_genome(s.assembly) if s.assembly else None
        log(f"    [{'OK     ' if fna else 'MISSING'}] genome FASTA {s.name}")
        if not fna:
            status["ok"] = False
            status["reasons"].append(f"no FASTA for {s.name}")

    if not status["ok"]:
        log("")
        log("    TIER 2 CANNOT RUN. These species will be EXCLUDED:", "ERROR")
        for s in needs:
            log(f"      - {s.name}"
                f"{'   <-- FOCAL SPECIES' if s.name == cfg.FOCAL_SPECIES else ''}",
                "ERROR")
        for r in status["reasons"]:
            log(f"      cause: {r}", "ERROR")
        log("    Fix: ln -sf $(conda info --base)/envs/acomys-bio/bin/"
            "{makeblastdb,blastn,blastdbcmd} $CONDA_PREFIX/bin/", "ERROR")
    else:
        log("    Tier 2 is ready.")
    log("")
    return status


def build(genes: list[str], species_list: list, dry_run: bool = False):
    cfg.ensure_dirs()
    tsa = cfg.REF / "transcriptome" / "acomys_15organ_TSA.fasta"
    tsa = tsa if tsa.exists() else None
    if tsa is None:
        log("15-organ TSA absent - A. cahirinus will fall back to PROJECTED "
            "annotation, which carries the circularity risk. Strongly "
            "recommend running 01_fetch_data.py --step transcriptome first.",
            "WARN")

    t2 = tier2_preflight(species_list) if not dry_run else {"ok": True}

    records: list[OrthologRecord] = []
    congener_cache: dict[str, str] = {}

    for gi, gene in enumerate(genes, 1):
        log(f"[{gi}/{len(genes)}] {gene}")
        if gene in cfg.ORTHOLOGY_WATCHLIST:
            log(f"    WATCHLIST: {cfg.ORTHOLOGY_WATCHLIST[gene]}", "WARN")

        # Resolve the congener first so its CDS can seed the TSA search.
        ordered = sorted(species_list,
                         key=lambda s: 0 if s.name == CONGENER_REFERENCE else 1)

        for sp in ordered:
            if dry_run:
                log(f"    DRY-RUN {sp.name}")
                continue
            rec = resolve_gene(gene, sp, tsa, congener_cache)
            records.append(rec)
            if sp.name == CONGENER_REFERENCE and rec.cds_length:
                fa = cfg.ORTHO / gene / f"{CONGENER_REFERENCE.replace(' ', '_')}.cds.fasta"
                if fa.exists():
                    congener_cache[gene] = str(next(SeqIO.parse(fa, "fasta")).seq)
            status = rec.confidence.upper()
            log(f"    {sp.name:24s} {rec.tier:14s} {rec.annotation_source:10s} "
                f"cds={rec.cds_length:5d} utr3={rec.utr3_length:5d} [{status}]")

    if dry_run or not records:
        return

    df = pd.DataFrame([asdict(r) for r in records])
    cfg.ORTHO.mkdir(parents=True, exist_ok=True)
    df.to_csv(cfg.ORTHO / "ortholog_table.csv", index=False)

    # ---- QC summary ----
    qc = (df.groupby("gene")
            .agg(n_species=("species", "nunique"),
                 n_high=("confidence", lambda s: (s == "high").sum()),
                 n_excluded=("confidence", lambda s: (s == "excluded").sum()),
                 n_projected=("annotation_source",
                              lambda s: (s == "projected").sum()),
                 median_cds=("cds_length", "median"),
                 median_utr3=("utr3_length", "median"))
            .reset_index())
    qc["usable_for_pgls"] = qc["n_species"] >= cfg.MIN_SPECIES_REQUIRED
    qc.to_csv(cfg.ORTHO / "qc_report.csv", index=False)

    log("\n" + "=" * 74)
    log("WP2 SUMMARY")
    log("=" * 74)
    log(f"gene x species records : {len(df)}")
    log(f"high confidence        : {(df.confidence == 'high').sum()}")
    log(f"medium                 : {(df.confidence == 'medium').sum()}")
    log(f"low                    : {(df.confidence == 'low').sum()}")
    log(f"excluded               : {(df.confidence == 'excluded').sum()}")

    log("")
    log("by tier   : " + ", ".join(f"{k}={v}" for k, v in
                                   df.tier.value_counts().items()))
    log("by source : " + ", ".join(f"{k}={v}" for k, v in
                                   df.annotation_source.value_counts().items()))

    # ---- the check that actually matters ----
    focal = df[df.species == cfg.FOCAL_SPECIES]
    n_focal_ok = (focal.confidence != "excluded").sum()
    log("")
    if n_focal_ok == 0:
        log(f"FOCAL SPECIES ({cfg.FOCAL_SPECIES}) HAS NO SEQUENCES.", "ERROR")
        log("WP3 measures Acomys deviation, so it CANNOT RUN. This is a hard "
            "stop, not a warning.", "ERROR")
        if not t2.get("ok", True):
            log("Cause: Tier 2 preflight failed - see the top of this run.",
                "ERROR")
        else:
            log("Tier 2 preflight passed, so homology recovery ran but found "
                "nothing. Check BLAST hits manually before proceeding.",
                "ERROR")
    else:
        log(f"focal species          : {n_focal_ok}/{len(focal)} genes "
            f"recovered", "INFO")

    excl = df[df.confidence == "excluded"]
    if len(excl):
        by_sp = excl.groupby("species").size().sort_values(ascending=False)
        log("")
        log("excluded by species:")
        for sp, n in by_sp.items():
            mark = "  <-- FOCAL" if sp == cfg.FOCAL_SPECIES else ""
            log(f"  {sp:24s} {n}{mark}")
        by_gene = excl.groupby("gene").size().sort_values(ascending=False)
        weak = by_gene[by_gene >= len(cfg.USABLE_SPECIES) - cfg.MIN_SPECIES_REQUIRED]
        if len(weak):
            log("")
            log("genes at risk of falling below the PGLS species minimum:")
            for g, n in weak.items():
                log(f"  {g:12s} missing in {n} species", "WARN")
    n_circ = (df.annotation_source == "projected").sum()
    if n_circ:
        log(f"\nPROJECTED annotation used in {n_circ} records - these carry the "
            f"circularity risk. Report WP3 with and without them.", "WARN")
    n_bad = (~qc.usable_for_pgls).sum()
    if n_bad:
        log(f"{n_bad} genes have < {cfg.MIN_SPECIES_REQUIRED} species and cannot "
            f"enter the PGLS.", "WARN")
    log(f"\nWrote: {cfg.ORTHO / 'ortholog_table.csv'}")
    log(f"       {cfg.ORTHO / 'qc_report.csv'}")
    log("Next : python 03_kmer_phylo_controlled.py --arm all")
    log("=" * 74)


def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--check", action="store_true",
                    help="report readiness and exit")
    ap.add_argument("--genes", choices=["focal", "controls", "all"],
                    default="focal")
    ap.add_argument("--all", action="store_true", help="same as --genes all")
    ap.add_argument("--dry-run", action="store_true")
    ap.add_argument("--skip-background", action="store_true")
    args = ap.parse_args()

    if args.check:
        sys.exit(0 if check_readiness() else 1)

    which_genes = "all" if args.all else args.genes
    if which_genes == "focal":
        genes = list(cfg.FOCAL_GENES)
    elif which_genes == "controls":
        genes = [g for s in cfg.CONTROL_SETS.values() for g in s]
    else:
        genes = cfg.all_genes()
        if not args.skip_background:
            genes += draw_random_background()

    log("=" * 74)
    log("WP2 - ORTHOLOG SCAFFOLD")
    log(f"genes   : {len(genes)}")
    log(f"species : {len(cfg.USABLE_SPECIES)}")
    log(f"pairs   : {len(genes) * len(cfg.USABLE_SPECIES)}")
    log("=" * 74)

    build(genes, cfg.USABLE_SPECIES, dry_run=args.dry_run)


if __name__ == "__main__":
    main()
