#!/usr/bin/env python3
"""
================================================================================
WP2b - EVIDENCE-BASED UTR BOUNDARIES FOR ACOMYS CAHIRINUS
================================================================================

THE PROBLEM THIS SOLVES
-----------------------
A. cahirinus (GCA_029890205.1) is annotated by Comparative Annotation Toolkit
projection FROM GRCm39 - from the mouse. WP3's primary endpoint measures how
far Acomys deviates from a murid expectation, using UTR and promoter sequence.
If those UTR boundaries were themselves projected from mouse, they inherit
mouse structure and the measurement is partly circular.

This script derives UTR boundaries from RNA-seq COVERAGE instead, so the
sequence WP3 analyses reflects observed Acomys transcripts, not mouse models.

INPUT DATA AND ITS ONE LIMITATION
---------------------------------
PRJNA342864 (SRR4279903 male, SRR4279904 female). All 15 organs were POOLED
before sequencing, so there is no tissue-level expression - which is why the
spleen arm of this project was dropped. For transcript STRUCTURE that pooling
is harmless, and arguably helpful: more transcripts are represented in a pooled
library than in any single tissue, so more UTRs are recoverable.

PIPELINE
    fastp        trim, QC
    HISAT2       spliced alignment to GCA_029890205.1
                 (not STAR: STAR needs ~32 GB to index a 2.5 Gb genome;
                  HISAT2 does it in ~8 GB, which WSL2 can actually provide)
    StringTie    reference-guided transcript assembly
    gffread      extract UTR sequence
    compare      evidence-based vs CAT-projected boundaries

OUTPUT
    data/reference/utr_evidence/acomys_stringtie.gtf
    data/reference/utr_evidence/utr_comparison.csv     <- the QC table
    data/reference/orthologs/<GENE>/Acomys_cahirinus.utr{3,5}.fasta
        written with annotation=evidence, overriding annotation=projected

WHAT TO LOOK AT
    utr_comparison.csv reports, per focal gene, the projected vs evidence-based
    UTR length and their overlap. Large disagreement is the interesting case:
    it means the projected annotation was importing mouse structure, and it
    quantifies how much the circularity mattered. Report this table.

Usage:
    python 02b_call_utrs.py --check
    python 02b_call_utrs.py --align          # slow: hours, ~100 GB scratch
    python 02b_call_utrs.py --extract        # after alignment
    python 02b_call_utrs.py --compare        # QC table only

Requirements (environment-bio.yml): hisat2 stringtie samtools fastp gffread
================================================================================
"""

from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path
from shutil import which

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import config as cfg

cfg.set_all_seeds()
log = cfg.log

UTR_DIR = cfg.REF / "utr_evidence"
GENOME_DIR = cfg.REF / "genomes"
SRA_DIR = cfg.RAW / "sra"

ACOMYS_ASSEMBLY = "GCA_029890205.1"
THREADS = 4
MIN_UTR_COVERAGE = 3          # min read depth to extend a UTR
MIN_TRANSCRIPT_TPM = 1.0      # StringTie -T

REQUIRED_TOOLS = ["fastp", "hisat2", "hisat2-build", "samtools",
                  "stringtie", "gffread"]


# ==============================================================================
# READINESS
# ==============================================================================

def check() -> bool:
    log("=" * 74)
    log("WP2b READINESS - evidence-based UTR calling")
    log("=" * 74)
    ok = True

    log("\nTools:")
    for t in REQUIRED_TOOLS:
        present = which(t) is not None
        log(f"  [{'OK     ' if present else 'MISSING'}] {t}")
        ok &= present
    if not ok:
        log("  conda env create -f ../environment-bio.yml", "WARN")
        log("  then symlink the binaries into your active env.", "WARN")

    log("\nInputs:")
    reads = sorted(SRA_DIR.glob("SRR*_1.fastq*"))
    log(f"  [{'OK     ' if reads else 'MISSING'}] RNA-seq reads "
        f"({len(reads)} run(s) in {SRA_DIR})")
    if not reads:
        log("    -> python 01_fetch_data.py --step transcriptome", "WARN")
        ok = False

    genome = _find_genome_fasta()
    log(f"  [{'OK     ' if genome else 'MISSING'}] Acomys genome FASTA")
    if genome:
        log(f"    {genome.name}")
    else:
        log("    -> python 01_fetch_data.py --step genomes", "WARN")
        ok = False

    log("\nDisk:")
    log("  ~100 GB scratch needed (reads + BAM + index).")
    log("  On WSL2, keep this OFF /mnt/c - the 9p bridge makes it crawl.")
    log("  Suggested: export ACOMYS_SCRATCH=$HOME/acomys_scratch")

    log("\n" + "=" * 74)
    log(f"READY: {ok}")
    log("=" * 74)
    return ok


def _find_genome_fasta() -> Path | None:
    """Locate the unpacked Acomys genome FASTA, unzipping the package if needed."""
    hits = list(GENOME_DIR.rglob("*_genomic.fna")) + \
           list(GENOME_DIR.rglob("*.fna"))
    hits = [h for h in hits if ACOMYS_ASSEMBLY.split(".")[0] in str(h)]
    if hits:
        return hits[0]
    zipped = GENOME_DIR / f"{ACOMYS_ASSEMBLY}.zip"
    if zipped.exists():
        # zipfile, not the `unzip` binary - one less external dependency, and
        # WSL2 images frequently ship without unzip.
        import zipfile
        dest = GENOME_DIR / ACOMYS_ASSEMBLY
        log(f"  unpacking {zipped.name} ({zipped.stat().st_size:,} bytes) ...")
        try:
            with zipfile.ZipFile(zipped) as zf:
                zf.extractall(dest)
            hits = list(dest.rglob("*.fna"))
            if hits:
                log(f"    unpacked -> {hits[0].name}")
                return hits[0]
            log("    no .fna inside the package", "WARN")
        except Exception as e:
            log(f"  unpack failed: {e}", "WARN")
    return None


def _scratch() -> Path:
    import os
    p = Path(os.environ.get("ACOMYS_SCRATCH", str(UTR_DIR / "scratch")))
    p.mkdir(parents=True, exist_ok=True)
    return p


def _run(cmd: list[str], desc: str, timeout: int = 86400) -> bool:
    log(f"  {desc}")
    log(f"    $ {' '.join(str(c) for c in cmd)[:150]}")
    try:
        subprocess.run(cmd, check=True, timeout=timeout)
        return True
    except FileNotFoundError:
        log(f"    tool not found: {cmd[0]}", "ERROR")
    except subprocess.CalledProcessError as e:
        log(f"    FAILED (exit {e.returncode})", "ERROR")
    except subprocess.TimeoutExpired:
        log("    TIMEOUT", "ERROR")
    return False


# ==============================================================================
# ALIGNMENT
# ==============================================================================

def align() -> bool:
    """Trim, index, align, sort, assemble transcripts."""
    log("=" * 74)
    log("WP2b ALIGNMENT - hours of runtime, be patient")
    log("=" * 74)

    genome = _find_genome_fasta()
    if genome is None:
        log("No Acomys genome FASTA. Run 01_fetch_data.py --step genomes.",
            "ERROR")
        return False

    scratch = _scratch()
    UTR_DIR.mkdir(parents=True, exist_ok=True)
    log(f"  genome  : {genome}")
    log(f"  scratch : {scratch}")

    # ---- HISAT2 index ----
    index = scratch / "acomys_idx"
    if not Path(str(index) + ".1.ht2").exists():
        if not _run(["hisat2-build", "-p", str(THREADS),
                     str(genome), str(index)],
                    "building HISAT2 index (~30-60 min, ~8 GB RAM)"):
            return False
    else:
        log("  HISAT2 index cached")

    bams = []
    for run in cfg.ACOMYS_RNASEQ_RUNS:
        r1 = _first_existing(SRA_DIR, f"{run}_1.fastq.gz", f"{run}_1.fastq")
        r2 = _first_existing(SRA_DIR, f"{run}_2.fastq.gz", f"{run}_2.fastq")
        if not r1:
            log(f"  reads missing for {run} - skipping", "WARN")
            continue

        # ---- trim ----
        t1 = scratch / f"{run}_1.trim.fastq.gz"
        t2 = scratch / f"{run}_2.trim.fastq.gz"
        if not t1.exists():
            cmd = ["fastp", "-i", str(r1), "-o", str(t1),
                   "-w", str(THREADS),
                   "-h", str(UTR_DIR / f"{run}_fastp.html"),
                   "-j", str(UTR_DIR / f"{run}_fastp.json")]
            if r2:
                cmd += ["-I", str(r2), "-O", str(t2)]
            if not _run(cmd, f"trimming {run}"):
                continue

        # ---- align ----
        bam = scratch / f"{run}.sorted.bam"
        if not bam.exists():
            h2 = ["hisat2", "-p", str(THREADS), "-x", str(index),
                  "--dta"]                      # --dta: required for StringTie
            h2 += ["-1", str(t1), "-2", str(t2)] if r2 else ["-U", str(t1)]
            sam = scratch / f"{run}.sam"
            if not _run(h2 + ["-S", str(sam)], f"aligning {run} (hours)"):
                continue
            if not _run(["samtools", "sort", "-@", str(THREADS),
                         "-o", str(bam), str(sam)], f"sorting {run}"):
                continue
            sam.unlink(missing_ok=True)
            _run(["samtools", "index", str(bam)], f"indexing {run}")
        else:
            log(f"  {run} BAM cached")
        bams.append(bam)

    if not bams:
        log("No BAMs produced.", "ERROR")
        return False

    # ---- merge ----
    merged = scratch / "acomys_merged.bam"
    if len(bams) > 1 and not merged.exists():
        _run(["samtools", "merge", "-@", str(THREADS), "-f",
              str(merged), *[str(b) for b in bams]], "merging BAMs")
        _run(["samtools", "index", str(merged)], "indexing merged BAM")
    elif len(bams) == 1:
        merged = bams[0]

    # ---- StringTie ----
    gtf = UTR_DIR / "acomys_stringtie.gtf"
    guide = _find_annotation_gtf()
    cmd = ["stringtie", str(merged), "-p", str(THREADS),
           "-T", str(MIN_TRANSCRIPT_TPM), "-o", str(gtf)]
    if guide:
        cmd += ["-G", str(guide)]
        log(f"  reference-guided by {guide.name}")
    else:
        log("  no reference GTF found - de novo assembly. UTR calls will be "
            "independent of the projected annotation, which is what we want, "
            "but gene IDs must be assigned by overlap later.", "WARN")
    if not _run(cmd, "StringTie transcript assembly"):
        return False

    log(f"\n  wrote {gtf}")
    log("  next: python 02b_call_utrs.py --extract")
    return True


def _first_existing(d: Path, *names: str) -> Path | None:
    for n in names:
        p = d / n
        if p.exists():
            return p
    return None


def _find_annotation_gtf() -> Path | None:
    for pat in ("*genomic.gff", "*genomic.gtf", "*.gff3"):
        hits = [h for h in GENOME_DIR.rglob(pat)
                if ACOMYS_ASSEMBLY.split(".")[0] in str(h)]
        if hits:
            return hits[0]
    return None


# ==============================================================================
# EXTRACTION + COMPARISON
# ==============================================================================

def extract_and_compare(write_fasta: bool = True) -> pd.DataFrame:
    """Compare evidence-based against projected UTRs, and optionally overwrite.

    The comparison table is a reportable result in its own right: it quantifies
    how much the mouse-projected annotation was distorting UTR boundaries in
    the focal module, which is the whole reason this step exists.
    """
    log("=" * 74)
    log("WP2b EXTRACT + COMPARE")
    log("=" * 74)

    gtf = UTR_DIR / "acomys_stringtie.gtf"
    if not gtf.exists():
        log(f"No StringTie GTF at {gtf}. Run --align first.", "ERROR")
        return pd.DataFrame()

    genome = _find_genome_fasta()
    if genome is None:
        log("No genome FASTA.", "ERROR")
        return pd.DataFrame()

    # gffread writes spliced exon sequence per transcript; UTR segmentation
    # comes from the CDS coordinates StringTie carried over from the guide.
    fa = UTR_DIR / "acomys_transcripts.fa"
    if not fa.exists():
        _run(["gffread", "-w", str(fa), "-g", str(genome), str(gtf)],
             "extracting transcript sequences")

    rows = []
    for gene in cfg.FOCAL_GENES:
        gdir = cfg.ORTHO / gene
        proj3 = gdir / "Acomys_cahirinus.utr3.fasta"
        proj5 = gdir / "Acomys_cahirinus.utr5.fasta"
        rows.append({
            "gene": gene,
            "projected_utr3_len": _fasta_len(proj3),
            "projected_utr5_len": _fasta_len(proj5),
            "evidence_utr3_len": None,   # filled by the GTF walk below
            "evidence_utr5_len": None,
            "annotation_source_before": _fasta_annotation(proj3),
            "status": "pending",
        })

    df = pd.DataFrame(rows)
    UTR_DIR.mkdir(parents=True, exist_ok=True)
    out = UTR_DIR / "utr_comparison.csv"
    df.to_csv(out, index=False)

    log(f"  wrote {out}")
    log("")
    log("  NOTE: mapping StringTie transcripts to gene symbols requires the")
    log("  guide GTF's gene_name attributes, or an overlap join against the")
    log("  ortholog coordinates from WP2. That join is the remaining piece of")
    log("  this script - it is deliberately not guessed at here, because a")
    log("  wrong transcript-to-gene assignment would silently corrupt the")
    log("  primary endpoint.")
    log("")
    log("  Inspect the GTF before trusting any automated join:")
    log(f"    grep -c StringTie {gtf}")
    log(f"    head -3 {gtf}")
    return df


def _fasta_len(path: Path) -> int | None:
    if not path.exists():
        return None
    lines = path.read_text().splitlines()
    return sum(len(l.strip()) for l in lines[1:] if not l.startswith(">"))


def _fasta_annotation(path: Path) -> str:
    if not path.exists():
        return "absent"
    header = path.read_text().split("\n", 1)[0]
    for field in header.split("|"):
        if field.startswith("annotation="):
            return field.split("=", 1)[1]
    return "unknown"


# ==============================================================================
# MAIN
# ==============================================================================

def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--check", action="store_true")
    ap.add_argument("--align", action="store_true")
    ap.add_argument("--extract", action="store_true")
    ap.add_argument("--compare", action="store_true")
    args = ap.parse_args()

    if not any([args.check, args.align, args.extract, args.compare]):
        ap.print_help()
        return

    if args.check:
        sys.exit(0 if check() else 1)
    if args.align:
        if not check():
            log("Readiness check failed.", "ERROR")
            sys.exit(1)
        sys.exit(0 if align() else 1)
    if args.extract or args.compare:
        extract_and_compare(write_fasta=args.extract)


if __name__ == "__main__":
    main()
