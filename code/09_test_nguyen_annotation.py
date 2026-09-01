#!/usr/bin/env python3
"""
================================================================================
DOES THE PUBLISHED A. cahirinus CAT ANNOTATION HAVE FRAME-BROKEN CDS?
================================================================================

This settles the single claim in Paper A that is currently UNSUPPORTED.

THE PROBLEM
-----------
We report that A. cahirinus coding sequence is frame-broken (internal stop
codons) in 27 of 35 genes, and we attribute this to CAT projection from GRCm39.

That attribution has never been tested. Our CDS did NOT come from a published
CAT annotation: NCBI publishes no annotation for GCA_029890205.1, so
01_fetch_data.py degraded to genome-only and 02_build_ortholog_scaffold.py
reconstructed the CDS by TIER-2 BLAST HSP STITCHING. Stitching can join HSPs
out of frame - we said so in 05_selection_analysis.py and then blamed CAT
anyway.

The observed pattern (A. russatus, native RefSeq, clean 34/35; A. cahirinus and
A. dimidiatus, tier-2, broken 27/35) is EQUALLY CONSISTENT with our own
stitching being the cause. Both explanations predict exactly what we saw.

THE TEST
--------
Nguyen ED et al. (2023) G3 13(10):jkad177 published a chromosome-scale
A. cahirinus assembly with a CAT annotation - the same method, done properly,
by the people who built the genome. Extract CDS for the same 35 genes from
their gff3 and run the identical reading-frame check.

  BROKEN TOO  -> our attribution is correct and the finding is STRONGER: it
                 concerns a published resource other groups are using, not an
                 artefact of our pipeline.
  CLEAN       -> the frame breaks are OURS. The CAT claim must be withdrawn
                 and replaced by a caution about tier-2 homology stitching.

Either result is publishable. Only the current untested state is not.

INPUTS (Zenodo 10.5281/zenodo.7761277, ~5.3 GB total)
------------------------------------------------------
  acahirinus.gff3               4.7 GB   CAT annotation
  acahirinus.scaffold.2bit      596 MB   genome sequence
Both are needed and must be used TOGETHER - scaffold names in the gff3 match
the 2bit, and will not necessarily match any NCBI assembly copy.

Usage:
    python code\\09_test_nguyen_annotation.py --download    # fetch from Zenodo
    python code\\09_test_nguyen_annotation.py --extract     # gff3 -> CDS
    python code\\09_test_nguyen_annotation.py --check       # frame test + verdict

Requires: py2bit  (pip install py2bit)
================================================================================
"""

from __future__ import annotations

import argparse
import gzip
import re
import struct
import sys
from collections import defaultdict
from pathlib import Path


# ==============================================================================
# MINIMAL 2bit READER  (pure Python - no compiler, no py2bit)
# ==============================================================================
# py2bit does not build on Windows: lib2bit/2bit.c includes <sys/mman.h>, which
# is POSIX-only, and there is no prebuilt wheel. The 2bit format is simple and
# stable, so it is read directly here rather than adding a toolchain dependency
# to a pipeline that already has to run natively on Windows.
#
# Format (UCSC): header {signature, version, seqCount, reserved} as uint32,
# then a name->offset index, then per sequence {dnaSize, nBlocks, maskBlocks,
# reserved, packedDna}. DNA is 2 bits per base, 4 bases per byte, most
# significant pair first, with T=0 C=1 A=2 G=3. Runs of N are stored separately
# as nBlocks and must be overlaid after unpacking - miss that step and every N
# silently reads as T, which would look like clean sequence and quietly corrupt
# a reading-frame test.

_BASES = "TCAG"
_SIG = 0x1A412743


class TwoBit:
    def __init__(self, path):
        self.fh = open(path, "rb")
        sig = struct.unpack("<I", self.fh.read(4))[0]
        if sig == _SIG:
            self.e = "<"
        else:
            self.fh.seek(0)
            sig = struct.unpack(">I", self.fh.read(4))[0]
            if sig != _SIG:
                raise ValueError(f"not a 2bit file (signature {sig:#x})")
            self.e = ">"
        _ver, n_seq, _res = struct.unpack(self.e + "III", self.fh.read(12))
        self.index = {}
        for _ in range(n_seq):
            ns = self.fh.read(1)[0]
            name = self.fh.read(ns).decode()
            off = struct.unpack(self.e + "I", self.fh.read(4))[0]
            self.index[name] = off
        self._hdr = {}

    def chroms(self):
        return dict.fromkeys(self.index)

    def _header(self, name):
        if name in self._hdr:
            return self._hdr[name]
        self.fh.seek(self.index[name])
        dna_size, n_cnt = struct.unpack(self.e + "II", self.fh.read(8))
        n_starts = struct.unpack(self.e + f"{n_cnt}I", self.fh.read(4 * n_cnt))
        n_sizes = struct.unpack(self.e + f"{n_cnt}I", self.fh.read(4 * n_cnt))
        m_cnt = struct.unpack(self.e + "I", self.fh.read(4))[0]
        self.fh.seek(8 * m_cnt + 4, 1)          # skip mask blocks + reserved
        h = (dna_size, list(zip(n_starts, n_sizes)), self.fh.tell())
        self._hdr[name] = h
        return h

    def sequence(self, name, start, end):
        dna_size, nblocks, dna_off = self._header(name)
        start = max(0, start); end = min(end, dna_size)
        if end <= start:
            return ""
        b0, b1 = start // 4, (end - 1) // 4
        self.fh.seek(dna_off + b0)
        raw = self.fh.read(b1 - b0 + 1)
        out = []
        for byte in raw:
            out.append(_BASES[(byte >> 6) & 3]); out.append(_BASES[(byte >> 4) & 3])
            out.append(_BASES[(byte >> 2) & 3]); out.append(_BASES[byte & 3])
        seq = out[start - b0 * 4: start - b0 * 4 + (end - start)]
        # Overlay N runs - these are NOT encoded in the 2-bit stream.
        for ns, nl in nblocks:
            ne = ns + nl
            if ne <= start or ns >= end:
                continue
            for i in range(max(ns, start), min(ne, end)):
                seq[i - start] = "N"
        return "".join(seq)

    def close(self):
        self.fh.close()

sys.path.insert(0, str(Path(__file__).resolve().parent))
import config as cfg

log = cfg.log
DEST = cfg.REF / "nguyen2023"
GFF = DEST / "acahirinus.gff3"
TWOBIT = DEST / "acahirinus.scaffold.2bit"
CDS_OUT = DEST / "nguyen_cds.fasta"
REPORT = cfg.RESULTS / "wp4_selection" / "nguyen_annotation_frame_check.csv"

URLS = {
    "acahirinus.gff3":
        "https://zenodo.org/records/7761277/files/acahirinus.gff3?download=1",
    "acahirinus.scaffold.2bit":
        "https://zenodo.org/records/7761277/files/acahirinus.scaffold.2bit?download=1",
}

# The same 35 genes tested in H2, so the comparison is like-for-like.
def target_genes() -> list[str]:
    g = {x.upper() for x in cfg.FOCAL_GENES}
    for v in cfg.CONTROL_SETS.values():
        g |= {x.upper() for x in (v or [])}
    return sorted(g)


def do_download():
    log("=" * 74)
    log("Download from Zenodo 10.5281/zenodo.7761277  (~5.3 GB)")
    log("=" * 74)
    DEST.mkdir(parents=True, exist_ok=True)
    log("\nThese are large. Run these two commands and let them finish:\n")
    for name, url in URLS.items():
        log(f'  curl -L -C - -o "{DEST / name}" "{url}"')
    log("\n  -C -  resumes a partial download rather than restarting.")
    log("\nThen:  python code\\09_test_nguyen_annotation.py --extract")
    return True


# ==============================================================================
# EXTRACT
# ==============================================================================

def _attr(s: str) -> dict:
    d = {}
    for kv in s.rstrip(";").split(";"):
        if "=" in kv:
            k, v = kv.split("=", 1)
            d[k.strip()] = v.strip()
    return d


def do_extract():
    """Stream the 4.7 GB gff3 once, keeping CDS features for target genes."""
    if not GFF.exists():
        log(f"missing {GFF} - run --download first", "ERROR")
        return False
    if not TWOBIT.exists():
        log(f"missing {TWOBIT} - run --download first", "ERROR")
        return False

    want = set(target_genes())
    # CAT writes the source symbol under varying keys depending on version;
    # try all of them rather than assuming one.
    NAME_KEYS = ("gene_name", "source_gene_common_name", "gene_id", "Name",
                 "gene", "Parent", "ID")

    # transcript -> list of (scaffold, start, end, strand, phase)
    tx_cds: dict[str, list] = defaultdict(list)
    tx_gene: dict[str, str] = {}
    seen_keys = defaultdict(int)

    opener = gzip.open if GFF.suffix == ".gz" else open
    log(f"streaming {GFF.name} ({GFF.stat().st_size/1e9:.1f} GB)...")
    n = 0
    with opener(GFF, "rt", errors="replace") as fh:
        for line in fh:
            n += 1
            if n % 5_000_000 == 0:
                log(f"  {n/1e6:.0f}M lines, {len(tx_cds)} transcripts kept")
            if not line or line[0] == "#":
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9 or f[2] != "CDS":
                continue
            # Cheap pre-filter before parsing attributes.
            up = f[8].upper()
            if not any(g in up for g in want):
                continue
            a = _attr(f[8])
            sym = None
            for k in NAME_KEYS:
                v = a.get(k)
                if v and v.upper() in want:
                    sym = v.upper(); seen_keys[k] += 1; break
            if sym is None:
                # symbol may be embedded, e.g. "Tgfb1-1"
                m = [g for g in want if re.search(rf"\b{g}\b", up)]
                if len(m) != 1:
                    continue
                sym = m[0]; seen_keys["regex"] += 1
            tid = a.get("transcript_id") or a.get("Parent") or a.get("ID")
            if not tid:
                continue
            tx_cds[tid].append((f[0], int(f[3]) - 1, int(f[4]), f[6],
                                f[7] if f[7] != "." else "0"))
            tx_gene[tid] = sym

    log(f"  done: {n:,} lines, {len(tx_cds)} transcripts across "
        f"{len(set(tx_gene.values()))} genes")
    log(f"  attribute keys used: {dict(seen_keys)}")

    tb = TwoBit(TWOBIT)
    chroms = set(tb.chroms())
    recs = {}
    for tid, parts in tx_cds.items():
        sym = tx_gene[tid]
        scaf = parts[0][0]
        if scaf not in chroms:
            continue
        strand = parts[0][3]
        parts = sorted(parts, key=lambda p: p[1], reverse=(strand == "-"))
        seq = []
        for s, st, en, _, _ in parts:
            try:
                sub = tb.sequence(s, st, en).upper()
            except Exception:
                sub = ""
            seq.append(sub)
        s = "".join(seq)
        if strand == "-":
            comp = str.maketrans("ACGTN", "TGCAN")
            s = s.translate(comp)[::-1]
        # keep the longest CDS per gene
        if sym not in recs or len(s) > len(recs[sym][1]):
            recs[sym] = (tid, s)
    tb.close()

    CDS_OUT.parent.mkdir(parents=True, exist_ok=True)
    with open(CDS_OUT, "w") as out:
        for sym, (tid, s) in sorted(recs.items()):
            out.write(f">{sym} {tid} len={len(s)}\n{s}\n")
    log(f"\nwrote {CDS_OUT}  ({len(recs)}/{len(want)} target genes recovered)")
    missing = sorted(want - set(recs))
    if missing:
        log(f"  not found in the annotation: {', '.join(missing)}", "WARN")
    return True


# ==============================================================================
# CHECK  - identical logic to 05_selection_analysis.prepare_cds
# ==============================================================================

def do_check():
    from Bio import SeqIO
    from Bio.Seq import Seq
    import pandas as pd

    if not CDS_OUT.exists():
        log(f"missing {CDS_OUT} - run --extract first", "ERROR")
        return False

    def frame_stats(nt: str):
        s = "".join(c if c in "ACGTacgt" else "N" for c in nt).upper()
        raw = len(s)
        t = len(s) % 3
        if t:
            s = s[:-t]
        pep = str(Seq(s).translate()) if s else ""
        if pep.endswith("*"):
            pep = pep[:-1]
        return raw, t == 0, pep.count("*")

    def scaffold_stops(gene: str, species: str):
        """Our own reconstruction, for the side-by-side comparison."""
        f = cfg.ORTHO / gene / f"{species}.cds.fasta"
        if not f.exists():
            return None
        try:
            return frame_stats(str(next(SeqIO.parse(f, "fasta")).seq))[2]
        except StopIteration:
            return None

    # All three tiers are written into ONE table, so downstream figures read a
    # single file instead of re-deriving the comparison from the scaffold. A
    # figure that recomputes its own inputs is a figure that can disagree with
    # the table it illustrates.
    rows = []
    for r in SeqIO.parse(CDS_OUT, "fasta"):
        raw_len, mult3, stops = frame_stats(str(r.seq))
        pep_amb = 0
        rows.append({
            "gene": r.id,
            "cds_len": raw_len,
            "multiple_of_3": mult3,
            "internal_stops": stops,
            "pct_ambiguous": pep_amb,
            "frame_clean": stops == 0 and mult3,
            "our_tier2_stops": scaffold_stops(r.id, "Acomys_cahirinus"),
            "russatus_refseq_stops": scaffold_stops(r.id, "Acomys_russatus"),
        })
    df = pd.DataFrame(rows).sort_values("internal_stops", ascending=False)
    REPORT.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(REPORT, index=False)

    log("=" * 74)
    log("FRAME CHECK — published CAT annotation (Nguyen et al. 2023)")
    log("=" * 74)
    log("\n" + df.to_string(index=False))

    n = len(df)
    broken = int((~df.frame_clean).sum())
    log("\n" + "=" * 74)
    log("THREE-TIER COMPARISON  (all computed here, none typed in)")
    log("=" * 74)
    for lab, col in (("native RefSeq (A. russatus)", "russatus_refseq_stops"),
                     ("published CAT (Nguyen 2023)", None),
                     ("homology recon. (this study)", "our_tier2_stops")):
        if col is None:
            b, t = broken, n
        else:
            sub = df[df[col].notna()]
            b, t = int((sub[col] > 0).sum()), len(sub)
        pct = 100 * b / t if t else float("nan")
        log(f"  {lab:30s} {b:2d}/{t:2d}  ({pct:.0f}%)")

    both = df[df.our_tier2_stops.notna()]
    if len(both):
        agree = ((both.internal_stops > 0) == (both.our_tier2_stops > 0)).mean()
        log(f"\n  concordance between the two reconstructions: {100*agree:.0f}% "
            f"({len(both)} shared genes)")
    log("=" * 74)

    if n == 0:
        log("  No genes recovered — check the attribute keys reported by "
            "--extract before drawing any conclusion.", "ERROR")
    elif broken / n > 0.4:
        log("\n  VERDICT: the published CAT annotation is ALSO frame-broken.")
        log("  Our attribution to CAT projection is SUPPORTED, and the finding")
        log("  is stronger than stated — it concerns a published resource in")
        log("  active use, not an artefact of our pipeline. Cite Nguyen et al.")
        log("  (2023) and report both figures side by side.")
    elif broken / n < 0.15:
        log("\n  VERDICT: the published CAT annotation is CLEAN.", "ERROR")
        log("  The frame breaks are OURS — tier-2 BLAST HSP stitching, not CAT", "ERROR")
        log("  projection. The CAT claim must be WITHDRAWN from both papers and", "ERROR")
        log("  replaced with a caution about homology-reconstructed CDS. H2", "ERROR")
        log("  should also be re-run using this annotation, which would restore", "ERROR")
        log("  a multi-species Deomyinae foreground.", "ERROR")
    else:
        log("\n  VERDICT: intermediate. Report the exact proportions and claim")
        log("  neither cause outright.", "WARN")
    log(f"\nwrote {REPORT}")
    return True


def main():
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--download", action="store_true")
    ap.add_argument("--extract", action="store_true")
    ap.add_argument("--check", action="store_true")
    a = ap.parse_args()
    if a.download:
        do_download()
    elif a.extract:
        sys.exit(0 if do_extract() else 1)
    elif a.check:
        sys.exit(0 if do_check() else 1)
    else:
        ap.print_help()


if __name__ == "__main__":
    main()
