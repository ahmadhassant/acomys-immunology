#!/usr/bin/env python3
"""
================================================================================
WP1 - DATA ACQUISITION
Acomys cahirinus fibrosis-free RENAL regeneration (kidney-only as of 13 Aug 2026)
================================================================================

Downloads and verifies every dataset the master plan depends on.

VERIFIED ACCESSIONS (checked 13 Aug 2026 against NCBI):
  GSE168876      Okamura & Majesky 2021, iScience. Acomys vs Mus kidney UUO.
                 18 samples = 2 species x {sham, UUO d2, UUO d5} x 3 reps.
                 BioProject PRJNA714406 | SRA SRP310563
                 Processed DE tables available as 4 supplementary XLSX files
                 -> no raw-read reprocessing needed for a first pass.
  GCA_029890205.1  ASM2989020v1. Nguyen et al. 2023 G3. Nanopore + Hi-C.
                   Scaffold N50 ~127 Mb, ~98.5% gene completeness. USE THIS.
  GCA_004027535.1  AcoCah_v1_BIUU. Superseded, fragmented. Fallback only.
  PRJNA342864      Mamrot et al. RNA-seq. WARNING: all 15 organs were POOLED
                   before sequencing (2 runs: SRR4279903 male, SRR4279904
                   female). NO tissue-level expression exists. Used ONLY for
                   evidence-based UTR boundaries. This is why the spleen arm
                   was dropped.
  GCF_030254825.1  Meriones unguiculatus Bangor_MerUng_6.1 (RefSeq annotated).
                   CRITICAL OUTGROUP: Gerbillinae is the sister subfamily to
                   Deomyinae (Acomys). Mus/Rattus are Murinae - more distant.

Usage:
    python 01_fetch_data.py --all
    python 01_fetch_data.py --step geo
    python 01_fetch_data.py --step orthologs --dry-run

Requirements:
    pip install biopython pandas requests openpyxl tqdm
    conda env create -f ../environment-bio.yml   # datasets, sra-tools, blast
    # EDirect is NOT required - Entrez work is done in Biopython
================================================================================
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import subprocess
import sys
import time
from pathlib import Path

import pandas as pd
import requests
from Bio import Entrez

# ==============================================================================
# CONFIGURATION - all constants live in config.py (single source of truth)
# ==============================================================================

sys.path.insert(0, str(Path(__file__).resolve().parent))
import config as cfg
import platform_compat as pc

cfg.set_all_seeds()

PROJECT_ROOT = cfg.PROJECT_ROOT
DATA, RAW, REF, LOGS = cfg.DATA, cfg.RAW, cfg.REF, cfg.LOGS

Entrez.email = cfg.ENTREZ_EMAIL
Entrez.api_key = cfg.ncbi_api_key()
REQUEST_DELAY = 0.12 if Entrez.api_key else 0.34

FOCAL_GENES = cfg.FOCAL_GENES
CONTROL_SETS = cfg.CONTROL_SETS
SPECIES_PANEL = cfg.SPECIES_PANEL
MIN_SPECIES_REQUIRED = cfg.MIN_SPECIES_REQUIRED
GSE168876_SAMPLES = cfg.GSE168876_SAMPLES

log = cfg.log
ensure_dirs = cfg.ensure_dirs

GEO_SUPPL_FILES = [
    "GSE168876_Acomys_day_2_vs_acomys_sham.xlsx",
    "GSE168876_Acomys_day_5_vs_acomys_sham.xlsx",
    "GSE168876_Mouse_day_2_vs_mouse_sham.xlsx",
    "GSE168876_Mouse_day_5_vs_mouse_sham.xlsx",
]
GEO_DL = "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE168876&format=file&file="


# ==============================================================================
# PROVENANCE MANIFEST
# ==============================================================================
# Every downloaded file is hashed and recorded. This is what makes the analysis
# reproducible and is worth its weight at revision time.

class Manifest:
    def __init__(self, path: Path):
        self.path = path
        self.entries: dict = {}
        if path.exists():
            self.entries = json.loads(path.read_text())

    def record(self, key: str, filepath: Path, source: str, note: str = ""):
        self.entries[key] = {
            "file": str(filepath.relative_to(PROJECT_ROOT)),
            "sha256": sha256(filepath),
            "bytes": filepath.stat().st_size,
            "source": source,
            "retrieved_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
            "note": note,
        }
        self.save()

    def save(self):
        self.path.parent.mkdir(parents=True, exist_ok=True)
        self.path.write_text(json.dumps(self.entries, indent=2))


def sha256(path: Path, chunk: int = 1 << 20) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        while block := fh.read(chunk):
            h.update(block)
    return h.hexdigest()


def have_tool(name: str) -> bool:
    from shutil import which
    return which(name) is not None


# ==============================================================================
# STEP 1 - GEO processed differential-expression tables
# ==============================================================================

def fetch_geo(manifest: Manifest, dry_run: bool = False):
    """Download the four processed DE tables. Fast path - no raw reads needed.

    NOTE the caveat that motivates step 3: the authors quantified Acomys reads
    against the MOUSE transcriptome. These tables inherit that bias. Use them
    for orientation and for reproducing the published result, then re-quantify
    natively before making any divergence claim.
    """
    log("STEP 1: GEO processed tables (GSE168876)")
    out = RAW / "geo"

    # Sample sheet
    sheet = out / "GSE168876_samples.csv"
    df = pd.DataFrame(GSE168876_SAMPLES,
                      columns=["gsm", "species", "condition", "day"])
    df["group"] = df["species"].str.split().str[0] + "_" + \
                  df["condition"] + "_d" + df["day"].astype(str)
    if not dry_run:
        df.to_csv(sheet, index=False)
        manifest.record("gse168876_samplesheet", sheet,
                        "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE168876",
                        "Transcribed from GEO record; 18 samples")
    log(f"  sample sheet -> {sheet} ({len(df)} samples)")

    for fname in GEO_SUPPL_FILES:
        dest = out / fname
        if dest.exists():
            log(f"  cached: {fname}")
            continue
        if dry_run:
            log(f"  DRY-RUN would fetch: {fname}")
            continue
        url = GEO_DL + fname.replace("_", "%5F").replace(".xlsx", "%2Exlsx")
        log(f"  fetching {fname} ...")
        try:
            r = requests.get(url, timeout=120)
            r.raise_for_status()
            dest.write_bytes(r.content)
            manifest.record(f"geo_{fname}", dest, url,
                            "Processed DE table (mouse-transcriptome mapped)")
            log(f"    ok ({len(r.content):,} bytes)")
        except Exception as e:
            log(f"    FAILED {fname}: {e}", "WARN")
        time.sleep(1.0)

    # Series matrix (metadata). GSE168876 spans TWO platforms, so GEO writes
    # one matrix per platform - there is no combined file.
    if not dry_run:
        base = ("https://ftp.ncbi.nlm.nih.gov/geo/series/GSE168nnn/"
                f"{cfg.GSE}/matrix/")
        got = []
        for gpl, desc in cfg.GSE_PLATFORMS.items():
            fname = f"{cfg.GSE}-{gpl}_series_matrix.txt.gz"
            dest = out / fname
            if dest.exists():
                log(f"  cached: {fname}")
                got.append(dest)
                continue
            try:
                r = requests.get(base + fname, timeout=120)
                r.raise_for_status()
                dest.write_bytes(r.content)
                manifest.record(f"geo_series_matrix_{gpl}", dest, base + fname,
                                f"Series metadata - {desc}")
                log(f"  series matrix {gpl} ok ({len(r.content):,} bytes)")
                got.append(dest)
            except Exception as e:
                log(f"  series matrix {gpl} failed: {e}", "WARN")
            time.sleep(0.5)
        if got:
            _verify_sample_sheet(got)


def _verify_sample_sheet(matrix_files: list[Path]):
    """Cross-check the hard-coded sample sheet against GEO's own metadata.

    cfg.GSE168876_SAMPLES was transcribed by hand from the GEO web record.
    Hand transcription is exactly the kind of thing that is silently wrong and
    poisons every downstream contrast, so verify it against the authoritative
    file rather than trusting it.
    """
    import gzip
    import re

    found: dict[str, str] = {}
    for path in matrix_files:
        try:
            with gzip.open(path, "rt", errors="replace") as fh:
                text = fh.read()
        except Exception as e:
            log(f"  could not read {path.name}: {e}", "WARN")
            continue
        gsms = re.findall(r"!Sample_geo_accession\s+(.+)", text)
        titles = re.findall(r"!Sample_title\s+(.+)", text)
        if not gsms:
            continue
        ids = [t.strip().strip('"') for t in gsms[0].split("\t") if t.strip()]
        names = ([t.strip().strip('"') for t in titles[0].split("\t") if t.strip()]
                 if titles else [""] * len(ids))
        for i, gsm in enumerate(ids):
            found[gsm] = names[i] if i < len(names) else ""

    expected = {g for g, _, _, _ in cfg.GSE168876_SAMPLES}
    actual = set(found)

    log("")
    log("  VERIFY sample sheet against GEO metadata:")
    log(f"    config: {len(expected)} samples | GEO: {len(actual)} samples")

    missing = expected - actual
    extra = actual - expected
    if not missing and not extra and len(expected) == len(actual):
        log("    MATCH - hand-transcribed sample sheet confirmed")
    else:
        if missing:
            log(f"    IN CONFIG BUT NOT IN GEO: {sorted(missing)}", "ERROR")
        if extra:
            log(f"    IN GEO BUT NOT IN CONFIG: {sorted(extra)}", "ERROR")
        log("    Fix cfg.GSE168876_SAMPLES before running WP5.", "ERROR")

    # Check species/condition/day labels too, not just IDs.
    mismatches = []
    for gsm, species, condition, day in cfg.GSE168876_SAMPLES:
        title = found.get(gsm, "").lower()
        if not title:
            continue
        sp_ok = species.split()[0].lower() in title
        cond_ok = (("sham" in title) if condition == "sham"
                   else ("uuo" in title))
        day_ok = (day == 0) or (f"day {day}" in title or f"day{day}" in title)
        if not (sp_ok and cond_ok and day_ok):
            mismatches.append((gsm, title))
    if mismatches:
        log(f"    LABEL MISMATCHES ({len(mismatches)}):", "ERROR")
        for gsm, title in mismatches[:6]:
            log(f"      {gsm}: GEO title = '{title}'", "ERROR")
    elif found:
        log("    Species / condition / day labels all consistent")


# ==============================================================================
# STEP 2 - Reference genomes
# ==============================================================================

def fetch_genomes(manifest: Manifest, dry_run: bool = False):
    """Pull assemblies via ncbi-datasets-cli, resolving accessions where unset."""
    log("STEP 2: reference genomes")
    if not pc.tool("datasets"):
        log("  'datasets' CLI not found. Install:", "WARN")
        log("    conda install -c conda-forge ncbi-datasets-cli", "WARN")
        return

    resolved: list[dict] = []
    for sp in SPECIES_PANEL:
        acc = sp.assembly or resolve_assembly(sp.taxid, sp.name)
        resolved.append({
            "species": sp.name, "taxid": sp.taxid, "clade": sp.clade,
            "regenerative": sp.regenerative, "assembly": acc,
            "assembly_name": sp.assembly_name, "annotation": sp.annotation,
            "usable": sp.usable, "note": sp.note,
        })
        if acc is None:
            log(f"  {sp.name}: NO ASSEMBLY FOUND - dropped", "WARN")
            continue
        if not sp.usable:
            log(f"  {sp.name}: {acc} - marked unusable ({sp.annotation}), "
                f"skipping download")
            continue
        log(f"  {sp.name}: {acc}  [{sp.clade}, {sp.annotation}]")
        if dry_run:
            continue
        dest = REF / "genomes" / f"{acc}.zip"
        if dest.exists() and dest.stat().st_size > 1000:
            log(f"    cached ({dest.stat().st_size:,} bytes)")
            continue

        # Unannotated assemblies (e.g. GCA_029890205.1) have no gff3/rna/cds/
        # protein, and `datasets` errors on the whole request rather than
        # skipping the missing parts. Degrade to genome-only.
        attempts = [("genome,gff3,rna,cds,protein", "full"),
                    ("genome", "genome-only")]
        for include, label in attempts:
            cmd = [pc.require("datasets"), "download", "genome", "accession", acc,
                   "--include", include, "--filename", str(dest)]
            try:
                subprocess.run(cmd, check=True, capture_output=True,
                               timeout=7200)
                manifest.record(f"genome_{acc}", dest,
                                f"NCBI Datasets {acc} ({include})",
                                f"{sp.name} - {sp.note}"[:200])
                log(f"    downloaded [{label}] "
                    f"({dest.stat().st_size:,} bytes)")
                if label == "genome-only":
                    log("    NOTE: no annotation in this package - WP2 must "
                        "derive gene models for this species.", "WARN")
                break
            except subprocess.CalledProcessError as e:
                msg = (e.stderr or b"").decode()[:200].strip()
                if label == "genome-only":
                    log(f"    FAILED: {msg}", "WARN")
                else:
                    log(f"    '{label}' failed, retrying genome-only "
                        f"({msg[:90]})")
            except subprocess.TimeoutExpired:
                log("    TIMEOUT after 2h", "WARN")
                break

    panel = REF / "genomes" / "species_panel_resolved.csv"
    if not dry_run:
        pd.DataFrame(resolved).to_csv(panel, index=False)

    n_ok = sum(1 for r in resolved if r["assembly"])
    log(f"  {n_ok}/{len(SPECIES_PANEL)} species have assemblies")
    if n_ok < MIN_SPECIES_REQUIRED:
        log(f"  Below the {MIN_SPECIES_REQUIRED}-species minimum for PGLS. "
            f"Fall back to the transcriptome-based ortholog route (step 4).", "WARN")


def resolve_assembly(taxid: int, name: str) -> str | None:
    """Find the best available assembly for a taxid. Prefers RefSeq (GCF)."""
    try:
        h = Entrez.esearch(db="assembly",
                           term=f"txid{taxid}[Organism:exp] AND "
                                f"(latest[filter] AND all[filter] "
                                f"NOT anomalous[filter])",
                           retmax=50)
        ids = Entrez.read(h)["IdList"]
        h.close()
        time.sleep(REQUEST_DELAY)
        if not ids:
            return None
        h = Entrez.esummary(db="assembly", id=",".join(ids))
        summ = Entrez.read(h)["DocumentSummarySet"]["DocumentSummary"]
        h.close()
        time.sleep(REQUEST_DELAY)

        def score(d):
            lvl = {"Chromosome": 3, "Complete Genome": 3,
                   "Scaffold": 2, "Contig": 1}.get(d.get("AssemblyStatus", ""), 0)
            refseq = 1 if d.get("RefSeq_category", "") != "na" else 0
            gcf = 1 if d.get("RefSeqAccession", "").startswith("GCF") else 0
            return (lvl, refseq, gcf)

        best = max(summ, key=score)
        return best.get("RefSeqAccession") or best.get("AssemblyAccession")
    except Exception as e:
        log(f"  resolve_assembly({name}): {e}", "WARN")
        return None


# ==============================================================================
# STEP 3 - Acomys 15-organ transcriptome (baseline spleen lives here)
# ==============================================================================

def fetch_transcriptome(manifest: Manifest, dry_run: bool = False):
    """PRJNA342864 raw reads - for UTR STRUCTURE, not expression.

    WHY RAW READS AND NOT A TSA ASSEMBLY
    ------------------------------------
    Mamrot et al. 2017 assembled 888,080 transcripts with EvidentialGene, but
    no TSA accession is given in the paper and none is findable in GenBank.
    Rather than chase it, we align the raw reads to the Acomys genome and call
    UTR boundaries from coverage. That is cheaper than re-deriving a de novo
    assembly and more direct for our purpose.

    WHAT THIS DATA IS AND IS NOT
    ----------------------------
    The 15 organs were POOLED before sequencing - two runs, male and female.
    So it carries NO tissue-level expression. That is fine here: we need
    transcript STRUCTURE (where UTRs start and end), and pooled tissue is
    arguably better for that, since more transcripts are represented.

    It is NOT a source of spleen expression. That is why the spleen arm was
    dropped (config.DESIGN_GAPS['no_spleen_data_at_all']).

    WHY IT MATTERS
    --------------
    A. cahirinus annotation is CAT-projected from GRCm39. Measuring divergence
    from mouse using mouse-projected UTR boundaries is circular. These reads
    give evidence-based boundaries instead.
    """
    log("STEP 3: Acomys RNA-seq for UTR evidence (PRJNA342864)")
    outdir = RAW / "sra"
    outdir.mkdir(parents=True, exist_ok=True)

    runs = cfg.ACOMYS_RNASEQ_RUNS
    log(f"  runs: {', '.join(runs)}  (pooled 15-organ; structure only)")

    if dry_run:
        for r in runs:
            log(f"  DRY-RUN would fetch {r} (~45 GB uncompressed each)")
        return

    if not pc.tool("fasterq-dump"):
        log("  sra-tools not found. Install:", "WARN")
        log("    conda env create -f ../environment-bio.yml", "WARN")
        log("  then symlink prefetch/fasterq-dump into your active env.", "WARN")
        _utr_manual_instructions()
        return

    # fasterq-dump writes tens of GB of temporary and output data. On WSL2 that
    # MUST go to Linux-native disk: writing across the /mnt/c 9p bridge fails
    # with "Cannot allocate memory(12)" partway through extraction.
    scratch = os.environ.get("ACOMYS_SCRATCH")
    if scratch:
        fq_out = Path(scratch) / "sra"
        tmpdir = Path(scratch) / "sra_tmp"
        fq_out.mkdir(parents=True, exist_ok=True)
        tmpdir.mkdir(parents=True, exist_ok=True)
        log(f"  extracting to Linux-native disk: {fq_out}")
    else:
        fq_out, tmpdir = outdir, outdir
        log("  ACOMYS_SCRATCH not set - extracting onto /mnt/c. On WSL2 this "
            "commonly fails with 'Cannot allocate memory'. Set ACOMYS_SCRATCH.",
            "WARN")

    import shutil as _shutil
    free_gb = _shutil.disk_usage(fq_out).free / 1e9
    log(f"  free space at target: {free_gb:.0f} GB")
    if free_gb < 200:
        log(f"  Extraction needs ~200 GB (.sra + FASTQ). Only {free_gb:.0f} GB "
            f"free. Consider --max-spots to subsample.", "WARN")

    for run in runs:
        r1 = fq_out / f"{run}_1.fastq.gz"
        if r1.exists() and r1.stat().st_size > 1_000_000:
            log(f"  cached: {run} ({r1.stat().st_size:,} bytes)")
            continue
        log(f"  fetching {run} - large (~35 GB .sra, ~90 GB FASTQ) and slow")
        try:
            subprocess.run(["prefetch", run, "--max-size", "100G",
                            "-O", str(outdir)],
                           check=True, timeout=36000)
            sra = _find_sra(outdir, run)
            # 2 threads, not 4: each thread holds its own buffers, and thread
            # count is the main lever on peak memory here.
            subprocess.run(["fasterq-dump", str(sra or run), "--split-files",
                            "--threads", "2", "--mem", "1000MB",
                            "-t", str(tmpdir), "-O", str(fq_out)],
                           check=True, timeout=36000)
            for suffix in ("_1", "_2"):
                fq = fq_out / f"{run}{suffix}.fastq"
                if fq.exists():
                    gz = pc.gzip_file(fq)   # no external gzip on Windows
                    if gz.exists():
                        manifest.record(f"sra_{run}{suffix}", gz,
                                        f"NCBI SRA {run}",
                                        cfg.ACOMYS_RNASEQ_NOTE)
            log(f"    {run} ok")
        except FileNotFoundError as e:
            log(f"    tool missing: {e}", "WARN")
            _utr_manual_instructions()
            return
        except subprocess.CalledProcessError as e:
            log(f"    FAILED {run}: {e}", "WARN")
        except subprocess.TimeoutExpired:
            log(f"    TIMEOUT on {run}", "WARN")


def _find_sra(outdir: Path, run: str) -> Path | None:
    """Locate the .sra prefetch already downloaded, so extraction reuses it."""
    for cand in (outdir / run / f"{run}.sra", outdir / f"{run}.sra"):
        if cand.exists():
            return cand
    hits = list(outdir.rglob(f"{run}.sra"))
    return hits[0] if hits else None


def _utr_manual_instructions():
    log("", "WARN")
    log("  UTR EVIDENCE IS NOT OPTIONAL. Without it, A. cahirinus falls back "
        "to mouse-PROJECTED annotation and WP3 measures divergence from mouse "
        "using boundaries derived from mouse - a circular result.", "WARN")
    log("  Manual route:", "WARN")
    log("    prefetch SRR4279903 SRR4279904", "WARN")
    log("    fasterq-dump --split-files SRR4279903 SRR4279904", "WARN")
    log("  Then: python 02b_call_utrs.py", "WARN")
    log("  Fallback if the reads are impractical: rely on A. russatus "
        "(GCF_903995435.1) RefSeq boundaries alone, and state the limitation.",
        "WARN")


# ==============================================================================
# STEP 4 - Ortholog sequences for the focal + control gene sets
# ==============================================================================

def fetch_orthologs(manifest: Manifest, dry_run: bool = False):
    """Pull CDS + mRNA (which carries the UTRs) per gene per species.

    UTRs matter more than CDS here: the regulatory hypothesis lives in 3'UTR
    AU-rich elements and miRNA seed sites, not in the strongly constrained
    coding sequence. We fetch the full mRNA so UTRs can be sliced out in WP2.
    """
    log("STEP 4: ortholog sequences")
    outdir = REF / "orthologs"
    all_genes = list(FOCAL_GENES) + [g for s in CONTROL_SETS.values() for g in s]
    log(f"  {len(all_genes)} genes x {len(SPECIES_PANEL)} species")

    index_rows = []
    for gene in all_genes:
        gdir = outdir / gene
        if not dry_run:
            gdir.mkdir(parents=True, exist_ok=True)
        for sp in SPECIES_PANEL:
            tag = sp.name.replace(" ", "_")
            fa = gdir / f"{tag}.mrna.fasta"
            if fa.exists():
                index_rows.append({"gene": gene, "species": sp.name,
                                   "file": str(fa), "status": "cached"})
                continue
            if dry_run:
                index_rows.append({"gene": gene, "species": sp.name,
                                   "file": str(fa), "status": "dry-run"})
                continue
            seqs = entrez_gene_mrna(gene, sp.taxid)
            if seqs:
                fa.write_text(seqs)
                index_rows.append({"gene": gene, "species": sp.name,
                                   "file": str(fa), "status": "ok"})
            else:
                index_rows.append({"gene": gene, "species": sp.name,
                                   "file": "", "status": "missing"})
            time.sleep(REQUEST_DELAY)
        log(f"  {gene} done")

    idx = outdir / "ortholog_index.csv"
    if not dry_run:
        pd.DataFrame(index_rows).to_csv(idx, index=False)
        manifest.record("ortholog_index", idx, "NCBI Gene/Nuccore",
                        f"{len(all_genes)} genes x {len(SPECIES_PANEL)} species")

    got = sum(1 for r in index_rows if r["status"] in ("ok", "cached"))
    log(f"  retrieved {got}/{len(index_rows)} gene-species pairs")
    log("  NOTE: gaps are expected for Acomys and the non-model species. "
        "02_build_ortholog_scaffold.py fills them by reciprocal-best-hit "
        "search against the assemblies and the de novo TSA.")


def entrez_gene_mrna(symbol: str, taxid: int) -> str | None:
    """Resolve gene symbol -> RefSeq mRNA (full length, UTRs included)."""
    try:
        h = Entrez.esearch(db="gene",
                           term=f"{symbol}[Gene Name] AND txid{taxid}[Organism]",
                           retmax=5)
        ids = Entrez.read(h)["IdList"]
        h.close()
        time.sleep(REQUEST_DELAY)
        if not ids:
            return None
        h = Entrez.elink(dbfrom="gene", db="nuccore", id=ids[0],
                         linkname="gene_nuccore_refseqrna")
        links = Entrez.read(h)
        h.close()
        time.sleep(REQUEST_DELAY)
        if not links or not links[0].get("LinkSetDb"):
            return None
        nuc = [l["Id"] for l in links[0]["LinkSetDb"][0]["Link"]][:3]
        h = Entrez.efetch(db="nuccore", id=",".join(nuc),
                          rettype="fasta", retmode="text")
        txt = h.read()
        h.close()
        return txt if txt.strip().startswith(">") else None
    except Exception:
        return None


# ==============================================================================
# STEP 5 - Cross-tissue replication datasets
# ==============================================================================

def list_replication_datasets():
    """Independent Acomys wound/regeneration datasets for WP6 replication.

    These are LEADS, not verified accessions - confirm each on GEO/SRA before
    committing. Cross-tissue replication is what turns a single-dataset result
    into a defensible one, given GSE168876 has only n=3 per group.
    """
    log("STEP 5: cross-tissue replication leads (verify manually)")
    leads = [
        ("Acomys dermal wound healing comparative transcriptome",
         "Brant/Maden et al., PLOS ONE 2019 - search GEO 'Acomys wound skin'"),
        ("Acomys skeletal muscle regeneration",
         "Maden et al., Sci Rep 2018"),
        ("Acomys vs Mus scar-free skin proteome",
         "Sci Rep 2019 - proteomic, use as orthogonal support"),
        ("Acomys single-cell muscle/dermis atlas",
         "Barbazuk R21-OD028209 - check for public release"),
        ("Mus musculus baseline spleen RNA-seq",
         "ENCODE / Mouse ENCODE - needed as the Mus arm of the spleen "
         "baseline comparison"),
    ]
    for name, note in leads:
        log(f"  - {name}\n      {note}")
    return leads


# ==============================================================================
# MAIN
# ==============================================================================

STEPS = {
    "geo": fetch_geo,
    "genomes": fetch_genomes,
    "transcriptome": fetch_transcriptome,
    "orthologs": fetch_orthologs,
}


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--all", action="store_true", help="run every step")
    ap.add_argument("--step", choices=list(STEPS) + ["replication"],
                    action="append", default=[])
    ap.add_argument("--dry-run", action="store_true",
                    help="report what would happen, download nothing")
    args = ap.parse_args()

    if not args.all and not args.step:
        ap.print_help()
        sys.exit(0)

    ensure_dirs()
    manifest = Manifest(DATA / "provenance_manifest.json")

    log("=" * 78)
    log("Acomys regeneration project - data acquisition")
    log(f"Project root : {PROJECT_ROOT}")
    log(f"NCBI API key : {'set' if Entrez.api_key else 'NOT SET (slower)'}")
    log(f"Mode         : {'DRY RUN' if args.dry_run else 'LIVE'}")
    log("=" * 78)

    todo = list(STEPS) if args.all else args.step
    for name in todo:
        if name == "replication":
            list_replication_datasets()
            continue
        try:
            STEPS[name](manifest, dry_run=args.dry_run)
        except Exception as e:
            log(f"step '{name}' raised: {e}", "ERROR")
        print()

    if args.all or "replication" in args.step:
        list_replication_datasets()

    log("=" * 78)
    log(f"Manifest: {manifest.path}")
    log(f"Entries : {len(manifest.entries)}")
    log("Next    : python 02_build_ortholog_scaffold.py")
    log("=" * 78)


if __name__ == "__main__":
    main()
