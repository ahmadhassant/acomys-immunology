#!/usr/bin/env python3
"""
================================================================================
WP5 - EXPRESSION RE-ANALYSIS (H3)
GSE168876 kidney UUO, Acomys cahirinus vs Mus musculus, sham / d2 / d5, n=3
================================================================================

H3 (PREREGISTRATION.md sec.3): Acomys constrains the TGF-beta1 and CCR2
pro-fibrotic arm after injury while preserving TGF-beta3 and CX3CL1-CX3CR1.

  Model            ~ species * timepoint;  contrast = species x time interaction
  PRIMARY          module-level test over curated fibrosis / macrophage sets
  SECONDARY        TGFB1:TGFB3 and CCR2:CX3CR1 ratio trajectories
  EXPLORATORY      gene-level DE - n=3 per group, labelled exploratory ALWAYS
  Falsified if     no significant species x time interaction at module level

WHY WE RE-QUANTIFY RATHER THAN USE THE PUBLISHED COUNTS
-------------------------------------------------------
The authors quantified Acomys reads against the MOUSE transcriptome.
Cross-species pseudo-alignment systematically under-recovers divergent genes,
and divergent genes are disproportionately chemokines and immune genes - i.e.
precisely our focal module. Using those counts would build the conclusion into
the input. This is carried design gap 'cross_species_mapping_bias'
(PREREGISTRATION.md sec.6 item 5); every divergence claim rests on species-native
re-quantification.

The direction of that bias matters: it would DEPLETE Acomys counts for the most
divergent focal genes, which could masquerade as "Acomys constrains this arm".
That is the H3 false positive we are specifically guarding against.

WHAT THIS ARM DOES *NOT* INHERIT
--------------------------------
  - H1's length confound       - this is expression, not sequence composition
  - H2's frame-break problem   - pseudo-alignment counts reads against
                                 transcripts; reading frame is irrelevant, so
                                 the CAT-projection defect does not propagate
  - H4's sort confound         - bulk tissue, no FACS gate
H3 is genuinely independent of every failure so far. That is the whole reason
it is still open.

POWER, STATED UP FRONT
----------------------
n = 3 per group. That is adequate for module-level tests aggregating hundreds
of genes and NOT adequate for gene-level DE. Gene-level output is labelled
exploratory in every file it appears in, and no gene-level claim will be made.

Usage:
    python code\\06_expression_reanalysis.py --check      # readiness, no compute
    python code\\06_expression_reanalysis.py --runs       # resolve GSM -> SRR
    python code\\06_expression_reanalysis.py --index      # build kallisto indices
    python code\\06_expression_reanalysis.py --quant      # quantify 18 runs
    python code\\06_expression_reanalysis.py --analyse    # H3 tests
================================================================================
"""

from __future__ import annotations

import argparse
import csv
import io
import json
import re
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import config as cfg
import platform_compat as pc

log = cfg.log

OUT = cfg.RESULTS / "wp5_expression"
QUANT = OUT / "quant"
IDX = cfg.REF / "kallisto"
FASTQ = cfg.RAW / "sra" / "GSE168876"
RUNS_JSON = cfg.RAW / "sra" / "GSE168876_runs.json"

# Species-native reference transcriptomes. Each species is quantified against
# ITS OWN transcriptome; nothing is mapped across species at the read level.
TRANSCRIPTOMES = {
    "Mus musculus": {
        "accession": "GCF_000001635.27",     # GRCm39
        "index": IDX / "Mus_musculus.idx",
        "fasta": cfg.REF / "transcriptome" / "Mus_musculus.rna.fna",
        "annotated": True,
    },
    # ------------------------------------------------------------------
    # A. cahirinus reads are quantified against the A. MINOUS transcriptome.
    #
    # WHY NOT A. cahirinus ITSELF: its assemblies (GCA_029890205.1,
    # GCA_030555615.1, GCA_004027535.1) are ALL unannotated GenBank entries.
    # NCBI ships no rna/cds/gff3 for any of them - 01_fetch_data.py already
    # degrades to genome-only for this accession. There is no A. cahirinus
    # transcript FASTA in existence to index.
    #
    # WHY A. minous IS THE RIGHT SUBSTITUTE, and not a compromise:
    # A. minous is not a distant congener - it sits INSIDE A. cahirinus
    # sensu lato. The Afro-Mediterranean clade (Cretan minous, Cypriot
    # nesiotes, Turkish cilicicus) is A. cahirinus in the broad sense;
    # karyotype banding and cytochrome b both place them very close, and
    # some authors treat them as varieties of one species. Splits within
    # the group are ~0.17 Ma (Giagia-Athanasopoulou et al. 2011, Biol J
    # Linn Soc 102:498; Barome et al., cahirinus-dimidiatus group).
    #
    # Scale of the improvement over the published counts:
    #     published (Acomys reads -> MOUSE)        ~18-20 Ma divergence
    #     A. russatus substitute                   ~5-10 Ma
    #     A. minous  (this)                        ~0.2-1 Ma
    # At that distance the residual mapping bias approaches within-species
    # polymorphism, which is the regime cross-species quantification needs
    # to be in before divergence claims are safe.
    #
    # This is still not literally "species-native" and IS a deviation from
    # the pre-registration. Declared in PREREGISTRATION.md sec.9.
    "Acomys cahirinus": {
        "accession": "GCF_964271855.1",      # A. minous mAcoMin1.4, RefSeq
        "reference_species": "Acomys minous",
        "index": IDX / "Acomys_minous.idx",
        "fasta": cfg.REF / "transcriptome" / "Acomys_minous.rna.fna",
        "annotated": True,
    },
}

# Candidate substitutes for the missing A. cahirinus transcriptome, best first.
ACOMYS_REFERENCE_OPTIONS = """
  1. A NEWER ANNOTATED A. cahirinus ASSEMBLY, if one now exists.
     Dissolves the problem entirely. Check first - it is free:
         python code\\check_new_data.py
  2. A. russatus RefSeq transcriptome (GCF_903995435.1).
     Congeneric substitute. Same genus, and the only Acomys with an
     independent RefSeq annotation - it is what rescued H2. Still
     CROSS-SPECIES, so it does not satisfy the pre-registration as written
     and needs a declared deviation. The residual bias runs in the same
     direction as the published mouse-based counts (divergent genes
     under-recovered) but is far smaller: ~2-10 Ma versus ~18-20 Ma.
  3. Genome-guided transcriptome assembly from the GSE168876 Acomys reads
     themselves (HISAT2/STAR + StringTie). Truly species-native and the
     methodologically correct answer, but needs a splice-aware aligner with
     no practical Windows build - this would have to run under WSL or Linux.
  4. Ortholog-scaffold index: quantify against the 241 recovered ortholog
     sequences only. Avoids cross-species mapping entirely, but a partial
     index breaks kallisto's EM assumption that the index spans the
     transcriptome - reads from the ~20,000 absent genes can pile onto
     included paralogs. NOT recommended as primary.
"""


# Samples excluded by --exclude, for sensitivity analyses. Empty by default;
# any exclusion is echoed in every downstream report so a filtered run can
# never be mistaken for the full one.
EXCLUDED_SAMPLES: set[str] = set()

# Canonical symbol -> symbol actually used in the Acomys minous annotation.
# Populated by --find-gene, which resolves genes that RefSeq left as LOC
# placeholders. Kept on disk so the mapping is auditable rather than hidden
# in code.
ALIAS_JSON = cfg.REF / "transcriptome" / "acomys_gene_aliases.json"


def acomys_aliases() -> dict[str, str]:
    if ALIAS_JSON.exists():
        return json.loads(ALIAS_JSON.read_text())
    return {}


def samples() -> list[tuple[str, str, str, int]]:
    return [s for s in cfg.GSE168876_SAMPLES if s[0] not in EXCLUDED_SAMPLES]


# ==============================================================================
# GSM -> SRR RESOLUTION
# ==============================================================================
# The 18 GSE168876 runs are NOT fetched by 01_fetch_data.py, which only handles
# PRJNA342864 (two runs, for UTR structure). Their SRR accessions are not in
# config because they are not printed in the GEO series record - each GSM links
# out to SRA individually.

RUNINFO_CSV = cfg.RAW / "sra" / "GSE168876_runinfo.csv"


def _entrez():
    from Bio import Entrez
    Entrez.email = cfg.ENTREZ_EMAIL
    Entrez.api_key = cfg.ncbi_api_key()
    return Entrez, (0.12 if Entrez.api_key else 0.34)


def _project_runinfo() -> str:
    """Whole-study runinfo CSV in one query.

    NOTE: 'GSM123[Accession]' does NOT work against db=sra. The Accession
    field there indexes SRA accessions (SRR/SRX/SRS/SRP) only - a GEO GSM is
    not in it, so the search returns zero hits for every sample. Querying the
    BioProject and mapping afterwards avoids the problem entirely, and costs
    one request instead of eighteen.
    """
    Entrez, delay = _entrez()
    for term in (f"{cfg.GSE_BIOPROJECT}[BioProject]",
                 cfg.GSE_BIOPROJECT, cfg.GSE_SRA):
        try:
            h = Entrez.esearch(db="sra", term=term, retmax=200)
            ids = Entrez.read(h)["IdList"]
            h.close()
            time.sleep(delay)
            if not ids:
                log(f"  '{term}' -> 0 records", "WARN")
                continue
            log(f"  '{term}' -> {len(ids)} SRA records")
            h = Entrez.efetch(db="sra", id=",".join(ids), rettype="runinfo",
                              retmode="text")
            txt = h.read()
            h.close()
            time.sleep(delay)
            if isinstance(txt, bytes):
                txt = txt.decode(errors="replace")
            if txt.strip():
                return txt
        except Exception as e:
            log(f"  '{term}' failed: {e}", "WARN")
    return ""


def _search_one(gsm: str) -> list[str]:
    """Fallback: free-text search. The GSM usually appears in SRA metadata as
    the sample alias, so a plain term hits even though [Accession] does not."""
    Entrez, delay = _entrez()
    try:
        h = Entrez.esearch(db="sra", term=gsm, retmax=20)
        ids = Entrez.read(h)["IdList"]
        h.close()
        time.sleep(delay)
        if not ids:
            return []
        h = Entrez.efetch(db="sra", id=",".join(ids), rettype="runinfo",
                          retmode="text")
        txt = h.read()
        h.close()
        time.sleep(delay)
        if isinstance(txt, bytes):
            txt = txt.decode(errors="replace")
        return sorted({r["Run"] for r in csv.DictReader(io.StringIO(txt))
                       if r.get("Run", "").startswith(("SRR", "ERR", "DRR"))})
    except Exception as e:
        log(f"  {gsm}: fallback search failed - {e}", "WARN")
        return []


def resolve_runs(force: bool = False) -> dict[str, list[str]]:
    """GSM -> [SRR, ...]. Cached to disk; NCBI is queried once."""
    if RUNS_JSON.exists() and not force:
        return json.loads(RUNS_JSON.read_text())

    RUNS_JSON.parent.mkdir(parents=True, exist_ok=True)
    out: dict[str, list[str]] = {}

    log("Strategy 1: whole-study runinfo via BioProject")
    txt = _project_runinfo()
    rows = []
    if txt:
        RUNINFO_CSV.write_text(txt)
        rows = [r for r in csv.DictReader(io.StringIO(txt)) if r.get("Run")]
        log(f"  {len(rows)} runs in runinfo -> {RUNINFO_CSV}")

    if rows:
        # Exact-token matching, NOT substring. 'GSM5171802' in 'GSM51718021'
        # is True, so a plain substring test can attach the wrong run to a
        # sample - which would silently mislabel its species or timepoint.
        def gsms_in(row) -> set[str]:
            found = set()
            for v in row.values():
                found.update(re.findall(r"GSM\d+", v or ""))
            return found

        for gsm, species, cond, day in samples():
            hits = sorted({r["Run"] for r in rows if gsm in gsms_in(r)})
            out[gsm] = hits
            log(f"  {gsm}  {species:17s} {cond:5s} d{day}  -> "
                f"{', '.join(hits) or 'unmatched'}")

        # No run may serve two samples.
        seen: dict[str, str] = {}
        for gsm, runs_ in out.items():
            for r in runs_:
                if r in seen:
                    log(f"  CONFLICT: {r} claimed by both {seen[r]} and "
                        f"{gsm} - refusing to guess", "ERROR")
                    out[gsm] = []
                    out[seen[r]] = []
                seen[r] = gsm

    # Only fall back for samples strategy 1 actually failed on. An earlier
    # version wrote `[...] or [all samples]`, and since an empty list is
    # falsy that re-queried all 18 even on a complete success.
    unmatched = ([s[0] for s in samples()] if not rows
                 else [g for g, v in out.items() if not v])
    if unmatched:
        log(f"\nStrategy 2: per-sample free-text search for "
            f"{len(unmatched)} unmatched", "WARN")
        for gsm in unmatched:
            out[gsm] = _search_one(gsm)
            log(f"  {gsm} -> {', '.join(out[gsm]) or 'NONE'}")

    n_ok = sum(1 for v in out.values() if v)
    RUNS_JSON.write_text(json.dumps(out, indent=2))
    log(f"\n{n_ok}/{len(samples())} samples resolved -> {RUNS_JSON}")

    if n_ok < len(samples()) and rows:
        # Make the failure diagnosable instead of just reporting "not found".
        log("\nUnresolved. What runinfo actually contains:", "WARN")
        log(f"  columns: {', '.join(list(rows[0].keys())[:18])}", "WARN")
        log("  first 3 rows, identifying fields only:", "WARN")
        keys = [k for k in ("Run", "Experiment", "SampleName", "BioSample",
                            "Sample", "LibraryName") if k in rows[0]]
        for r in rows[:3]:
            log("    " + "  ".join(f"{k}={r.get(k,'')}" for k in keys), "WARN")
        log(f"\n  Full CSV is at {RUNINFO_CSV} - paste the identifying "
            "columns and the GSM mapping can be made explicit.", "WARN")
    elif n_ok < len(samples()):
        log("\nNo runinfo retrieved at all. Check the NCBI API key and "
            "network, then retry with --runs --force.", "WARN")
    return out


def fastq_for(srr: str) -> list[Path]:
    """Reads on disk for one run, paired or single, gzipped or not."""
    hits = []
    for pat in (f"{srr}_1.fastq*", f"{srr}_2.fastq*", f"{srr}.fastq*"):
        hits.extend(sorted(FASTQ.glob(pat)))
    return hits


# ==============================================================================
# READINESS
# ==============================================================================

def check() -> bool:
    log("=" * 74)
    log("WP5 H3 READINESS - expression re-analysis")
    log("=" * 74)
    ok = True

    log("\n1. DESIGN")
    by = {}
    for gsm, sp, cond, day in samples():
        by.setdefault((sp, cond, day), []).append(gsm)
    for (sp, cond, day), g in sorted(by.items()):
        log(f"   {sp:17s} {cond:5s} d{day}  n={len(g)}")
    log(f"   total {len(samples())} samples, 2 species x 3 conditions x n=3")
    log("   n=3 -> module-level tests only. Gene-level DE is exploratory.")

    log("\n2. TOOLS")
    for t in ("kallisto", "prefetch", "fasterq-dump"):
        p = pc.tool(t)
        log(f"   [{'OK     ' if p else 'MISSING'}] {t}  {p or ''}")
        if t == "kallisto" and not p:
            ok = False
    if not pc.tool("kallisto"):
        # The GitHub releases page no longer publishes Windows assets. The
        # official SourceForge mirror still does - use it rather than
        # concluding kallisto is Linux-only.
        log("   kallisto (Windows) - NOT on the GitHub releases page any "
            "more; use the official mirror:", "WARN")
        log("     https://sourceforge.net/projects/kallisto.mirror/files/"
            "v0.51.1/kallisto_windows-v0.51.1.tar.gz/download", "WARN")
        log("     older fallback: .../files/v0.48.0/kallisto_windows-"
            "v0.48.0.zip/download", "WARN")
        log("   extract to <drive>:\\Tools\\kallisto, or set ACOMYS_TOOL_DIRS",
            "WARN")
    if not pc.tool("fasterq-dump"):
        log("   sra-toolkit (Windows): sratoolkit.<ver>-win64.zip from", "WARN")
        log("     https://github.com/ncbi/sra-tools/wiki/"
            "01.-Downloading-SRA-Toolkit", "WARN")
        log("   extract to <drive>:\\Tools\\sratoolkit-<ver>", "WARN")

    log("\n3. REFERENCE TRANSCRIPTOMES  (species-native - the point of this arm)")
    blocked = False
    for sp, meta in TRANSCRIPTOMES.items():
        fa, ix = meta["fasta"], meta["index"]
        has_fa = fa.exists() and fa.stat().st_size > 1_000_000
        has_ix = ix.exists()
        size = f"{fa.stat().st_size/1e6:.0f} MB" if fa.exists() else "-"
        state = "OK" if has_fa else ("NOT OBTAINABLE" if not meta["annotated"]
                                     else "missing, downloadable")
        via = meta.get("reference_species")
        tag = f"{sp} (via {via})" if via else sp
        log(f"   {tag:38s} {meta['accession']:18s} "
            f"fasta {state:22s} {size:>8s}  "
            f"index {'OK' if has_ix else 'missing'}")
        if not has_fa:
            ok = False
            if not meta["annotated"]:
                blocked = True

    if blocked:
        log("\n   BLOCKER - read this before downloading anything.", "ERROR")
        log("   GCA_029890205.1 is an UNANNOTATED GenBank assembly. NCBI "
            "publishes no transcript FASTA for it, so there is nothing to "
            "build an A. cahirinus kallisto index from.", "ERROR")
        log("   'Species-native re-quantification' as pre-registered is "
            "therefore not directly executable. Substitutes:", "ERROR")
        for line in ACOMYS_REFERENCE_OPTIONS.strip("\n").splitlines():
            log(f"   {line}", "WARN")
        log("\n   Do NOT start the ~150 GB read download until this is "
            "settled - the choice determines whether those reads are "
            "usable.", "ERROR")

    log("\n4. READS")
    runs = json.loads(RUNS_JSON.read_text()) if RUNS_JSON.exists() else {}
    if not runs:
        log("   SRR accessions not resolved yet -> run --runs", "WARN")
        ok = False
    else:
        n_have = n_missing = 0
        for gsm, sp, cond, day in samples():
            srrs = runs.get(gsm, [])
            present = [s for s in srrs if fastq_for(s)]
            if srrs and len(present) == len(srrs):
                n_have += 1
            else:
                n_missing += 1
                log(f"   {gsm} {sp:17s} {cond:5s} d{day}  "
                    f"missing {[s for s in srrs if s not in present] or 'SRR unknown'}",
                    "WARN")
        log(f"   {n_have}/{len(samples())} samples have reads on disk")
        if n_missing:
            ok = False

    log("\n5. ORTHOLOG SPACE")
    ortho = cfg.ORTHO
    n_genes = len(list(ortho.glob("*"))) if ortho.exists() else 0
    log(f"   scaffold: {n_genes} genes at {ortho}")
    if n_genes < 100:
        log("   run 02_build_ortholog_scaffold.py --all first", "WARN")
        ok = False

    log("\n" + "=" * 74)
    log(f"READY: {ok}")
    if not ok:
        log("\nOrder of operations:")
        log("  1. --runs    resolve GSM -> SRR (fast, one NCBI query)")
        log("  2. download reads with prefetch + fasterq-dump")
        log("  3. fetch both transcriptomes")
        log("  4. --index   build kallisto indices")
        log("  5. --quant   quantify all 18 runs species-natively")
        log("  6. --analyse H3 tests")
    log("=" * 74)
    return ok


def fetch_reference():
    """Download both RefSeq transcriptomes via ncbi-datasets, unzip the
    rna.fna. Both are GCF_ (RefSeq) accessions, so an annotation - and
    therefore rna.fna - is expected; the code verifies rather than assumes."""
    import subprocess
    import zipfile

    datasets = pc.tool("datasets")
    if datasets is None:
        log("'datasets' CLI not found.", "ERROR")
        log("  https://www.ncbi.nlm.nih.gov/datasets/docs/v2/download-and-install/",
            "WARN")
        return False

    dest_dir = cfg.REF / "transcriptome"
    dest_dir.mkdir(parents=True, exist_ok=True)
    ok = True

    for sp, meta in TRANSCRIPTOMES.items():
        fa = meta["fasta"]
        if fa.exists() and fa.stat().st_size > 1_000_000:
            log(f"  {sp}: cached ({fa.stat().st_size/1e6:.0f} MB)")
            continue
        acc = meta["accession"]
        zp = dest_dir / f"{acc}.zip"
        log(f"  {sp}: downloading {acc} (rna)")
        try:
            subprocess.run([datasets, "download", "genome", "accession", acc,
                            "--include", "rna", "--filename", str(zp)],
                           check=True, capture_output=True, timeout=7200)
        except Exception as e:
            err = getattr(e, "stderr", b"")
            err = err.decode(errors="replace")[:300] if err else str(e)
            log(f"    FAILED: {err}", "ERROR")
            log("    If this reports no RNA in the package, the assembly is "
                "not annotated and the reference choice must be revisited.",
                "ERROR")
            ok = False
            continue

        with zipfile.ZipFile(zp) as z:
            rna = [n for n in z.namelist() if n.endswith("rna.fna")]
            if not rna:
                log(f"    no rna.fna inside {zp.name} - assembly carries no "
                    "transcripts", "ERROR")
                ok = False
                continue
            with z.open(rna[0]) as src, open(fa, "wb") as out:
                out.write(src.read())
        n_tx = sum(1 for line in fa.read_text(errors="replace").splitlines()
                   if line.startswith(">"))
        log(f"    -> {fa.name}  {fa.stat().st_size/1e6:.0f} MB, "
            f"{n_tx:,} transcripts")
        if n_tx < 10_000:
            log(f"    WARNING: only {n_tx:,} transcripts. A mammalian "
                "transcriptome should have tens of thousands; a sparse "
                "index silently depresses mapping rate.", "WARN")
    return ok


# ==============================================================================
# INDEX + QUANTIFY
# ==============================================================================
# The 18 runs are roughly 100-180 GB of FASTQ if extracted all at once. They are
# NOT all extracted at once. Each run is fetched, quantified, and deleted before
# the next begins, so peak disk stays around 15 GB.
#
# This is not premature tidiness. An earlier fasterq-dump run in this project
# defaulted its temp directory to the current working directory on a mounted
# Windows drive and orphaned 101 GB. Scratch location and per-run cleanup are
# therefore explicit here.

def read_run_info(path: Path) -> dict:
    """Parse kallisto's run_info.json, tolerating malformed Windows output.

    kallisto writes the invoking command line verbatim into the JSON 'call'
    field. On Windows that contains backslashes - D:\\Tools\\kallisto\\
    kallisto.exe - which are not escaped, so '\\T' and '\\k' are invalid JSON
    escapes and a strict parser rejects the whole file. The quantification is
    fine; only the metadata is malformed. Falling over here would discard a
    completed run, so parse defensively.
    """
    txt = path.read_text(errors="replace")
    try:
        return json.loads(txt)
    except json.JSONDecodeError:
        pass
    # Escape any backslash that is not already a valid JSON escape.
    try:
        return json.loads(re.sub(r'\\(?!["\\/bfnrtu])', r"\\\\", txt))
    except json.JSONDecodeError:
        pass
    # Last resort: pull the numbers out directly.
    out: dict = {}
    for key in ("n_targets", "n_bootstraps", "n_processed",
                "n_pseudoaligned", "n_unique", "p_pseudoaligned",
                "p_unique"):
        m = re.search(rf'"{key}"\s*:\s*"?(-?[\d.]+)"?', txt)
        if m:
            v = m.group(1)
            out[key] = float(v) if "." in v else int(v)
    if not out:
        log(f"  could not parse {path}", "WARN")
    return out


def build_indices() -> bool:
    import subprocess
    kal = pc.tool("kallisto")
    if kal is None:
        log("kallisto not found", "ERROR")
        return False
    IDX.mkdir(parents=True, exist_ok=True)
    for sp, meta in TRANSCRIPTOMES.items():
        ix, fa = meta["index"], meta["fasta"]
        if ix.exists():
            log(f"  {sp}: index cached")
            continue
        if not fa.exists():
            log(f"  {sp}: {fa} missing - run --fetch-ref", "ERROR")
            return False
        log(f"  {sp}: indexing {fa.name} (needs ~8-16 GB RAM, several minutes)")
        try:
            pc.run_tool(kal, ["index", "-i", str(ix), str(fa)],
                        capture_output=True, text=True, timeout=14400)
            log(f"    -> {ix.name}  {ix.stat().st_size/1e6:.0f} MB")
        except Exception as e:
            out = getattr(e, "stderr", "") or getattr(e, "stdout", "") or str(e)
            log(f"    FAILED: {str(out)[-400:]}", "ERROR")
            return False
    return True


def quant_all(keep_fastq: bool = False) -> bool:
    """Fetch -> quantify -> delete, one run at a time."""
    import shutil
    import subprocess

    kal = pc.tool("kallisto")
    fqd = pc.tool("fasterq-dump")
    pre = pc.tool("prefetch")
    if not (kal and fqd):
        log("kallisto and fasterq-dump are both required", "ERROR")
        return False

    runs = resolve_runs()
    QUANT.mkdir(parents=True, exist_ok=True)
    scratch = pc.space_free_dir(cfg.RAW / "sra" / "scratch")
    rows = []

    for gsm, species, cond, day in samples():
        srrs = runs.get(gsm, [])
        if not srrs:
            log(f"  {gsm}: no SRR - skipped", "WARN")
            continue
        srr = srrs[0]
        odir = QUANT / f"{gsm}_{srr}"
        info = odir / "run_info.json"

        if info.exists():
            log(f"  {gsm} {srr}: cached")
        else:
            odir.mkdir(parents=True, exist_ok=True)
            FASTQ.mkdir(parents=True, exist_ok=True)
            try:
                if pre:
                    log(f"  {gsm} {srr}: prefetch")
                    subprocess.run([pre, srr, "--max-size", "100G",
                                    "--output-directory", str(FASTQ)],
                                   check=True, capture_output=True, timeout=14400)
                log(f"  {gsm} {srr}: extracting FASTQ")
                subprocess.run([fqd, srr, "--split-files", "--skip-technical",
                                "-e", "2", "--mem", "1000MB",
                                "-t", str(scratch), "-O", str(FASTQ)],
                               check=True, capture_output=True, timeout=14400,
                               cwd=str(FASTQ))
            except Exception as e:
                err = getattr(e, "stderr", b"")
                err = err.decode(errors="replace")[:300] if err else str(e)
                log(f"    FAILED to fetch: {err}", "ERROR")
                continue

            fq = fastq_for(srr)
            if not fq:
                log(f"    no FASTQ produced for {srr}", "ERROR")
                continue

            idx = TRANSCRIPTOMES[species]["index"]
            paired = len([p for p in fq if p.name.endswith(("_1.fastq",
                                                            "_1.fastq.gz",
                                                            "_2.fastq",
                                                            "_2.fastq.gz"))]) >= 2
            args = ["quant", "-i", str(idx), "-o", str(odir), "-t", "4"]
            if paired:
                args += [str(fq[0]), str(fq[1])]
            else:
                # Fragment length must be supplied for single-end; these are
                # nominal values and are recorded so the choice is auditable.
                args += ["--single", "-l", "200", "-s", "30", str(fq[0])]
            log(f"    quantifying ({'paired' if paired else 'single'}-end) "
                f"vs {idx.name}")
            try:
                pc.run_tool(kal, args, capture_output=True, text=True,
                            timeout=14400)
            except Exception as e:
                out = getattr(e, "stderr", "") or str(e)
                log(f"    kallisto FAILED: {str(out)[-400:]}", "ERROR")
                continue
            finally:
                if not keep_fastq:
                    for p in fastq_for(srr):
                        p.unlink(missing_ok=True)
                    sra = FASTQ / srr
                    if sra.is_dir():
                        shutil.rmtree(sra, ignore_errors=True)

        if info.exists():
            ri = read_run_info(info)
            rows.append({
                "gsm": gsm, "srr": srr, "species": species,
                "condition": cond, "day": day,
                "n_processed": ri.get("n_processed"),
                "n_pseudoaligned": ri.get("n_pseudoaligned"),
                "p_pseudoaligned": ri.get("p_pseudoaligned"),
            })

    if not rows:
        log("no samples quantified", "ERROR")
        return False

    import pandas as pd
    df = pd.DataFrame(rows)
    df.to_csv(OUT / "mapping_rates.csv", index=False)

    log("\n" + "=" * 74)
    log("MAPPING RATE BY SPECIES - the audit of the A. minous substitution")
    log("=" * 74)
    log("\n" + df.groupby("species").p_pseudoaligned
        .agg(["count", "mean", "min", "max"]).round(2).to_string())

    means = df.groupby("species").p_pseudoaligned.mean()
    if len(means) == 2:
        gap = abs(means.iloc[0] - means.iloc[1])
        log(f"\nspecies gap: {gap:.1f} percentage points")
        if gap > 15:
            log("WARNING: large asymmetry. A lower Acomys mapping rate means "
                "residual reference bias survives, and it depletes exactly "
                "the divergent genes the focal module is made of - the same "
                "direction as the published mouse-based counts, just "
                "smaller. Report this alongside every H3 result and treat "
                "any 'Acomys constrains this arm' finding as suspect.",
                "ERROR")
        else:
            log("Comparable across species - the substitution is not "
                "introducing an obvious asymmetric bias.")
    log(f"\nwrote {OUT / 'mapping_rates.csv'}")
    return True


# ==============================================================================
# H3 ANALYSIS
# ==============================================================================

def tx2gene(fasta: Path) -> dict[str, str]:
    """RefSeq transcript accession -> gene symbol, from FASTA headers.

    RefSeq deflines look like
        >NM_007393.5 Mus musculus actin, beta (Actb), mRNA
    The symbol is the LAST parenthesised group. Taking the last one rather
    than the first matters for names that themselves contain brackets, e.g.
    'solute carrier family 2 (facilitated glucose transporter), member 1
    (Slc2a1)' - the first group is prose, the last is the symbol.

    SYMBOLS ARE UPPERCASED. The two references use different case
    conventions - mouse RefSeq is Titlecase (Actb, Zfp85), the Acomys minous
    annotation is ALLCAPS (ACTB, CASKIN2; 28,630 of 28,631 symbols). Matching
    on raw case gave an overlap of FIFTEEN genes out of ~29,000 and ~40,000.
    Uppercasing both sides gives 15,447. Config gene sets are already
    uppercase; the marker panels imported from 04 are Titlecase, so they are
    uppercased at the point of use too.
    """
    out: dict[str, str] = {}
    with open(fasta, "r", errors="replace") as fh:
        for line in fh:
            if not line.startswith(">"):
                continue
            head = line[1:].rstrip()
            acc, _, desc = head.partition(" ")
            groups = re.findall(r"\(([^()]+)\)", desc)
            if groups:
                out[acc] = groups[-1].strip().upper()
    return out


def diagnose_symbols(n: int = 12):
    """Why do the two symbol namespaces not intersect?

    A shared-symbol count in the tens rather than the thousands means the two
    references are not speaking the same naming language. The usual causes,
    which this distinguishes:
      - one annotation uses LOC######## placeholders for uncurated genes
      - case conventions differ (ACTB / Actb / actb)
      - the defline format differs, so the parser is picking up prose
    """
    log("=" * 74)
    log("SYMBOL NAMESPACE DIAGNOSTIC")
    log("=" * 74)
    parsed = {}
    for sp, meta in TRANSCRIPTOMES.items():
        fa = meta["fasta"]
        if not fa.exists():
            log(f"{sp}: {fa} missing", "ERROR")
            continue
        label = meta.get("reference_species", sp)
        log(f"\n--- {label}  ({fa.name}) ---")
        heads = []
        with open(fa, "r", errors="replace") as fh:
            for line in fh:
                if line.startswith(">"):
                    heads.append(line[1:].rstrip())
                    if len(heads) >= n:
                        break
        for h in heads:
            acc, _, desc = h.partition(" ")
            groups = re.findall(r"\(([^()]+)\)", desc)
            sym = groups[-1].strip() if groups else None
            log(f"  {acc:18s} sym={str(sym):18s} | {desc[:70]}")

        t2g = tx2gene(fa)
        syms = set(t2g.values())
        parsed[sp] = syms
        loc = sum(1 for s in syms if s.upper().startswith("LOC"))
        upper = sum(1 for s in syms if s.isupper())
        lower = sum(1 for s in syms if s.islower())
        title = sum(1 for s in syms if s[:1].isupper() and not s.isupper())
        log(f"\n  transcripts with a symbol : {len(t2g):,}")
        log(f"  distinct symbols          : {len(syms):,}")
        log(f"  LOC-prefixed              : {loc:,} ({100*loc/max(len(syms),1):.1f}%)")
        log(f"  ALLCAPS / Titlecase / lower: {upper:,} / {title:,} / {lower:,}")

    if len(parsed) == 2:
        a, b = list(parsed.values())
        log("\n" + "=" * 74)
        log("OVERLAP")
        log("=" * 74)
        log(f"  exact                      : {len(a & b):,}")
        ua = {s.upper() for s in a}
        ub = {s.upper() for s in b}
        log(f"  case-insensitive           : {len(ua & ub):,}")
        na = {s for s in a if not s.upper().startswith("LOC")}
        nb = {s for s in b if not s.upper().startswith("LOC")}
        log(f"  case-insensitive, no LOC   : "
            f"{len({s.upper() for s in na} & {s.upper() for s in nb}):,}")
        for want in ("TGFB1", "TGFB3", "CCR2", "CX3CR1", "CX3CL1", "ACTB"):
            log(f"    {want:8s} A={want in ua}  B={want in ub}")
    return True


def find_gene(symbol: str, top: int = 8) -> bool:
    """Locate a gene in the Acomys minous transcriptome by homology.

    Genes that NCBI's pipeline could not confidently assign to an ortholog are
    named LOC########, so they never match by symbol - 19.8% of A. minous
    symbols are LOC-prefixed. CCR2 is one of them. Rapidly evolving families
    with lineage-specific duplications are exactly the ones this happens to,
    which is also why CCL2 was dropped at WP2 orthology QC.

    Strategy: take the mouse transcript(s) for the symbol, blastn against the
    Acomys transcriptome, and report the best hits with their symbols. The
    mapping is NOT applied automatically - a reciprocal-best-hit judgement on
    a chemokine receptor family is exactly where a paralog gets mistaken for
    an ortholog, so the decision is left explicit.
    """
    import subprocess
    import tempfile

    blastn = pc.tool("blastn")
    makedb = pc.tool("makeblastdb")
    if not (blastn and makedb):
        log("blastn/makeblastdb not found", "ERROR")
        return False

    symbol = symbol.upper()
    mus_fa = TRANSCRIPTOMES["Mus musculus"]["fasta"]
    aco_fa = TRANSCRIPTOMES["Acomys cahirinus"]["fasta"]

    # Pull the mouse transcripts for this symbol.
    t2g_mus = tx2gene(mus_fa)
    want = {acc for acc, s in t2g_mus.items() if s == symbol}
    if not want:
        log(f"{symbol} not found in the mouse transcriptome", "ERROR")
        return False
    log(f"{symbol}: {len(want)} mouse transcript(s)")

    tmp = Path(tempfile.mkdtemp())
    query = tmp / f"{symbol}.fna"
    keep, buf, cur = [], [], None
    with open(mus_fa, "r", errors="replace") as fh:
        for line in fh:
            if line.startswith(">"):
                if cur in want and buf:
                    keep.append("".join(buf))
                cur = line[1:].split()[0]
                buf = [line] if cur in want else []
            elif buf:
                buf.append(line)
    if cur in want and buf:
        keep.append("".join(buf))
    query.write_text("".join(keep))

    db = aco_fa.with_suffix(".blastdb")
    if not Path(str(db) + ".nin").exists() and not Path(str(db) + ".nal").exists():
        log("  building BLAST db for the Acomys transcriptome (once)")
        subprocess.run([makedb, "-in", pc.safe_path(aco_fa), "-dbtype", "nucl",
                        "-out", pc.safe_path(db)],
                       check=True, capture_output=True, timeout=7200)

    log("  blastn ...")
    r = subprocess.run(
        [blastn, "-query", pc.safe_path(query), "-db", pc.safe_path(db),
         "-outfmt", "6 qseqid sseqid pident length evalue bitscore",
         "-max_target_seqs", "20", "-evalue", "1e-20"],
        check=True, capture_output=True, text=True, timeout=7200)

    t2g_aco = tx2gene(aco_fa)
    seen, rows = set(), []
    for line in r.stdout.splitlines():
        f = line.split("\t")
        if len(f) < 6:
            continue
        sym = t2g_aco.get(f[1], "?")
        if sym in seen:
            continue
        seen.add(sym)
        rows.append((f[1], sym, float(f[2]), int(f[3]), f[4], float(f[5])))
        if len(rows) >= top:
            break

    if not rows:
        log(f"  no hits for {symbol} - it may genuinely be absent", "WARN")
        return False

    log(f"\n  {'subject':20s} {'symbol':18s} {'pident':>7s} {'len':>6s} "
        f"{'evalue':>10s} {'bits':>8s}")
    for acc, sym, pid, ln, ev, bs in rows:
        log(f"  {acc:20s} {sym:18s} {pid:7.1f} {ln:6d} {ev:>10s} {bs:8.1f}")

    best = rows[0]
    log(f"\n  Best forward hit: {best[1]}  ({best[2]:.1f}% identity over "
        f"{best[3]} bp, bits {best[5]:.0f})")
    if len(rows) > 1:
        ratio = best[5] / rows[1][5]
        log(f"  Margin over runner-up ({rows[1][1]}): {ratio:.2f}x bitscore")

    # ------------------------------------------------------------------
    # RECIPROCAL BEST HIT - the test that actually establishes orthology
    # ------------------------------------------------------------------
    # A forward best hit is NOT evidence of orthology. Blast the candidate
    # back against the mouse transcriptome: if it returns the gene we started
    # from, the two are reciprocal best hits. If it returns a neighbour, the
    # forward hit was the closest available paralog, not the ortholog.
    #
    # This matters acutely here. CCR2 sits in a chemokine receptor cluster
    # with CCR1/CCR3/CCR5, and CCR5 is CCR2's immediate neighbour - it turned
    # up as the second forward hit, with HIGHER percent identity than the
    # winner. Accepting a forward hit in this family is precisely the mistake
    # that put a wrong ortholog into WP2 and got CCL2 removed.
    log("\n  RECIPROCAL CHECK - blasting the candidate back against mouse")
    cand_acc = best[0]
    cand_fa = tmp / "candidate.fna"
    buf, cur, got = [], None, False
    with open(aco_fa, "r", errors="replace") as fh:
        for line in fh:
            if line.startswith(">"):
                if got:
                    break
                cur = line[1:].split()[0]
                if cur == cand_acc:
                    got, buf = True, [line]
            elif got:
                buf.append(line)
    cand_fa.write_text("".join(buf))

    mdb = mus_fa.with_suffix(".blastdb")
    if not (Path(str(mdb) + ".nin").exists()
            or Path(str(mdb) + ".nal").exists()):
        log("    building BLAST db for the mouse transcriptome (once)")
        subprocess.run([makedb, "-in", pc.safe_path(mus_fa), "-dbtype", "nucl",
                        "-out", pc.safe_path(mdb)],
                       check=True, capture_output=True, timeout=7200)

    rb = subprocess.run(
        [blastn, "-query", pc.safe_path(cand_fa), "-db", pc.safe_path(mdb),
         "-outfmt", "6 sseqid pident length evalue bitscore",
         "-max_target_seqs", "20", "-evalue", "1e-20"],
        check=True, capture_output=True, text=True, timeout=7200)

    back, seen_b = [], set()
    for line in rb.stdout.splitlines():
        f = line.split("\t")
        if len(f) < 5:
            continue
        s = t2g_mus.get(f[0], "?")
        if s in seen_b:
            continue
        seen_b.add(s)
        back.append((s, float(f[1]), int(f[2]), float(f[4])))
        if len(back) >= 6:
            break

    log(f"\n  {'mouse symbol':18s} {'pident':>7s} {'len':>6s} {'bits':>8s}")
    for s, pid, ln, bs in back:
        log(f"  {s:18s} {pid:7.1f} {ln:6d} {bs:8.1f}")

    ok = bool(back) and back[0][0] == symbol
    log("")
    if ok:
        margin = back[0][3] / back[1][3] if len(back) > 1 else float("inf")
        log(f"  RECIPROCAL BEST HIT CONFIRMED: {cand_acc} <-> {symbol} "
            f"(margin {margin:.2f}x over {back[1][0] if len(back)>1 else 'nothing'})")
        log(f"  Safe to add:  {{\"{symbol}\": \"{best[1]}\"}}  ->  {ALIAS_JSON}")
    else:
        top = back[0][0] if back else "nothing"
        log(f"  RECIPROCAL CHECK FAILED. {cand_acc} maps back to {top}, "
            f"not {symbol}.", "ERROR")
        log(f"  The forward hit is the nearest available PARALOG, not the "
            f"{symbol} ortholog. Do NOT add the alias - a wrong ortholog is "
            "more damaging than a missing one, because it produces a "
            "plausible-looking expression profile.", "ERROR")
    return ok


def gene_tpm(species: str) -> "pd.DataFrame":
    """genes x samples TPM for one species, summed over transcripts."""
    import pandas as pd
    runs = resolve_runs()
    t2g = tx2gene(TRANSCRIPTOMES[species]["fasta"])
    cols = {}
    for gsm, sp, cond, day in samples():
        if sp != species:
            continue
        srrs = runs.get(gsm, [])
        if not srrs:
            continue
        ab = QUANT / f"{gsm}_{srrs[0]}" / "abundance.tsv"
        if not ab.exists():
            log(f"  missing {ab}", "WARN")
            continue
        df = pd.read_csv(ab, sep="\t")
        df["gene"] = df["target_id"].map(t2g)
        df = df.dropna(subset=["gene"])
        cols[gsm] = df.groupby("gene")["tpm"].sum()
    if not cols:
        return pd.DataFrame()
    out = pd.DataFrame(cols)

    # Rename LOC placeholders to their canonical symbols where --find-gene has
    # established the correspondence by homology.
    if species == "Acomys cahirinus":
        alias = acomys_aliases()
        rev = {v: k for k, v in alias.items()}
        hits = [g for g in out.index if g in rev]
        if hits:
            out = out.rename(index=rev)
            log(f"  applied {len(hits)} Acomys alias(es): "
                f"{', '.join(f'{h}->{rev[h]}' for h in hits)}")
    return out


def modules() -> dict[str, list[str]]:
    """Gene sets tested at module level. All symbols UPPERCASED to match
    tx2gene() - config sets are already uppercase, the marker panels are not."""
    mods = {"FOCAL": [g.upper() for g in cfg.FOCAL_GENES]}

    # Secondary endpoint: E. Hassanat's a priori panel, three blocks tested
    # SEPARATELY so the boundary of the effect is located rather than averaged
    # away, with the combined 15 and full 20 reported alongside.
    for k, v in cfg.EXTENDED_PANEL.items():
        mods[f"PANEL_{k}"] = [g.upper() for g in v]
    mods["PANEL_combined15"] = [g.upper() for g in cfg.EXTENDED_15]
    mods["PANEL_full20"] = [g.upper() for g in cfg.PANEL_20]

    for k, v in cfg.CONTROL_SETS.items():
        if v:
            mods[k] = [g.upper() for g in v]
    # Macrophage-polarisation panels, as pre-registered. Reused verbatim from
    # the single-cell arm so the two arms cannot drift apart.
    try:
        import importlib.util
        p = Path(__file__).with_name("04_scrna_myeloid.py")
        spec = importlib.util.spec_from_file_location("sc", p)
        sc = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(sc)
        for k in ("macro_proinflam", "macro_resolving", "fibrogenic",
                  "antifibrotic", "monocyte_classical", "monocyte_patrol"):
            if k in sc.MARKERS:
                mods[k] = [g.upper() for g in sc.MARKERS[k]]
    except Exception as e:
        log(f"  could not import marker panels from 04: {e}", "WARN")
    return mods


def analyse(n_perm: int = 10000, seed: int = 42) -> bool:
    import numpy as np
    import pandas as pd

    OUT.mkdir(parents=True, exist_ok=True)
    log("=" * 74)
    log("H3 - EXPRESSION DIVERGENCE (species x time interaction)")
    log("=" * 74)

    aco = gene_tpm("Acomys cahirinus")
    mus = gene_tpm("Mus musculus")
    if aco.empty or mus.empty:
        log("quantification missing - run --quant first", "ERROR")
        return False

    # INNER join on shared symbols. An outer join inflated a cross-species
    # comparison earlier in this project (WP4-sc); genes absent from one
    # reference must not be scored as zero in that species.
    shared = sorted(set(aco.index) & set(mus.index))
    log(f"\nGenes: Acomys {len(aco.index):,}  Mus {len(mus.index):,}  "
        f"SHARED {len(shared):,}")
    if len(shared) < 5000:
        log("Shared-symbol count is implausibly low for two mammalian "
            "annotations. Run --symbols before trusting anything below.",
            "ERROR")
    aco, mus = aco.loc[shared], mus.loc[shared]

    # Which focal genes actually survived into shared space, and where each
    # missing one was lost. A module quietly shrinking from 5 genes to 4
    # changes what the primary endpoint means.
    log("\nFocal module membership in shared space:")
    shared_set = set(shared)
    dropped = []
    for g in (x.upper() for x in cfg.FOCAL_GENES):
        where = ("shared" if g in shared_set else
                 "Mus only" if g in mus.index or g in set(mus.index) else
                 "absent from both")
        if g not in shared_set:
            dropped.append(g)
        log(f"  {g:10s} {where}")
    if dropped:
        log(f"\n  {len(dropped)} focal gene(s) lost: {', '.join(dropped)}",
            "WARN")
        log("  CCR2 in particular is absent from the A. minous annotation. "
            "The CCL2-CCR2 chemokine cluster is rapidly evolving with "
            "lineage-specific duplications, which is the same reason CCL2 "
            "was dropped from the module at WP2 orthology QC. It is likely "
            "present in the assembly under a LOC placeholder rather than "
            "truly absent - 19.8% of A. minous symbols are LOC-prefixed.",
            "WARN")
        log("  Consequence: the focal module is tested with "
            f"{len(cfg.FOCAL_GENES) - len(dropped)} genes, not "
            f"{len(cfg.FOCAL_GENES)}. State this with the result.", "WARN")

    meta = {g: (sp, cond, day) for g, sp, cond, day in samples()}

    def cols_for(df, cond, day):
        return [c for c in df.columns
                if meta[c][1] == cond and meta[c][2] == day]

    # log2FC vs that species' OWN sham. Each species is its own baseline, so
    # constant per-species reference bias cancels in the contrast.
    def lfc(df, day):
        sham = df[cols_for(df, "sham", 0)].mean(axis=1)
        inj = df[cols_for(df, "UUO", day)].mean(axis=1)
        return np.log2(inj + 1) - np.log2(sham + 1)

    res = {}
    for day in (2, 5):
        a, m = lfc(aco, day), lfc(mus, day)
        res[day] = pd.DataFrame({"lfc_acomys": a, "lfc_mus": m,
                                 "interaction": a - m})

    rng = np.random.default_rng(seed)
    mods = modules()
    rows = []
    for day, df in res.items():
        stat_all = df["interaction"].values
        for name, genes in mods.items():
            present = [g for g in genes if g in df.index]
            if len(present) < 4:
                log(f"  {name} d{day}: only {len(present)} genes present "
                    "- skipped", "WARN")
                continue
            obs = float(df.loc[present, "interaction"].mean())
            # COMPETITIVE test: is this module's interaction shifted relative
            # to all other genes? Null = random gene sets of equal size, which
            # preserves the observed interaction distribution.
            k = len(present)
            null = np.array([rng.choice(stat_all, k, replace=False).mean()
                             for _ in range(n_perm)])
            p_two = float((np.abs(null - null.mean())
                           >= abs(obs - null.mean())).mean())
            rows.append({"day": day, "module": name, "n_genes": k,
                         "mean_interaction": obs,
                         "null_mean": float(null.mean()),
                         "null_sd": float(null.std()),
                         "z": float((obs - null.mean()) / (null.std() or np.nan)),
                         "p_perm": p_two})

    tab = pd.DataFrame(rows)
    if tab.empty:
        log("no module had enough genes", "ERROR")
        return False
    tab["q_BH"] = _bh(tab["p_perm"].values)
    tab = tab.sort_values(["day", "p_perm"])
    tab.to_csv(OUT / "h3_module_tests.csv", index=False)

    # Secondary endpoint: report which panel genes resolved. A gene that fails
    # to resolve is REPORTED, never substituted - the panel is fixed.
    log("\nEXTENDED PANEL resolution in shared space (secondary endpoint):")
    for blk, genes in cfg.EXTENDED_PANEL.items():
        present = [g for g in genes if g in shared_set]
        missing = [g for g in genes if g not in shared_set]
        log(f"  {blk:24s} {len(present)}/{len(genes)}"
            + (f"   MISSING: {', '.join(missing)}" if missing else ""))
    log("\nPRIMARY - module-level species x time interaction")
    log("(interaction = Acomys log2FC - Mus log2FC, vs each species' own sham)")
    log("\n" + tab.to_string(index=False))

    log("\n" + "-" * 74)
    focal = tab[tab.module == "FOCAL"]
    sig = focal[focal.q_BH < cfg.ALPHA]
    log("PRE-SPECIFIED RULE: H3 is falsified if there is no significant "
        "species x time interaction at module level.")
    if len(sig):
        log(f"RESULT: FOCAL module significant at {len(sig)} timepoint(s) "
            f"-> H3 SUPPORTED", )
    else:
        log("RESULT: FOCAL module not significant at any timepoint "
            "-> H3 NOT SUPPORTED")

    # ------------------------------------------------------------------
    # MANDATORY CONFOUND CHECK - the direct analogue of H1's length covariate
    # ------------------------------------------------------------------
    # The observed pattern (divergent immune genes shifted, conserved
    # housekeeping genes null) is what the hypothesis predicts. It is ALSO
    # exactly what relative read recovery would produce on its own:
    #
    #   A gene recovering fewer reads in Acomys has its log2FC compressed
    #   toward zero by the +1 pseudocount floor. Compressed Acomys log2FC
    #   minus intact Mus log2FC = NEGATIVE interaction. Conserved genes map
    #   equally well in both species and show no such compression.
    #
    # So recovery imbalance alone predicts: focal negative, housekeeping null.
    # H1 produced a beautiful region-specific 3'UTR signal that turned out to
    # be sequence length. This is the same test, applied before the claim.
    #
    # The model is fitted on NON-FOCAL genes only, so focal membership cannot
    # influence the correction - the same discipline used in WP3.
    from scipy import stats as sstats

    log("\n" + "=" * 74)
    log("CONFOUND CHECK - is the interaction explained by read recovery?")
    log("=" * 74)

    a_mean = np.log2(aco.mean(axis=1) + 1)
    m_mean = np.log2(mus.mean(axis=1) + 1)
    recovery = a_mean - m_mean          # <0 => gene recovers less in Acomys
    abundance = (a_mean + m_mean) / 2

    focal_set = {g.upper() for g in cfg.FOCAL_GENES}
    conf_rows = []
    for day, df in res.items():
        inter = df["interaction"]
        nonfocal = [g for g in df.index if g not in focal_set]

        rho_r, p_r = sstats.spearmanr(recovery.loc[nonfocal],
                                      inter.loc[nonfocal])
        rho_a, p_a = sstats.spearmanr(abundance.loc[nonfocal],
                                      inter.loc[nonfocal])
        log(f"\nday {day}, within NON-FOCAL genes (n={len(nonfocal):,}):")
        log(f"  interaction vs relative recovery : Spearman {rho_r:+.3f}  "
            f"p={p_r:.3g}")
        log(f"  interaction vs mean abundance    : Spearman {rho_a:+.3f}  "
            f"p={p_a:.3g}")

        fg = [g for g in focal_set if g in df.index]
        log(f"  focal genes' recovery percentile : "
            f"{[round(float(sstats.percentileofscore(recovery.loc[nonfocal], recovery.loc[g])), 1) for g in fg]}")
        log(f"    (50 = typical; low values mean focal genes recover poorly "
            f"in Acomys relative to Mus)")

        # Residualise interaction on recovery + abundance, MODEL FITTED ON
        # NON-FOCAL GENES ONLY, then re-test every module.
        X = np.column_stack([np.ones(len(nonfocal)),
                             recovery.loc[nonfocal].values,
                             abundance.loc[nonfocal].values])
        y = inter.loc[nonfocal].values
        beta, *_ = np.linalg.lstsq(X, y, rcond=None)
        Xall = np.column_stack([np.ones(len(df)), recovery.loc[df.index].values,
                                abundance.loc[df.index].values])
        resid = pd.Series(inter.values - Xall @ beta, index=df.index)

        stat_all = resid.values
        for name, genes in mods.items():
            present = [g for g in genes if g in resid.index]
            if len(present) < 4:
                continue
            obs = float(resid.loc[present].mean())
            k = len(present)
            null = np.array([rng.choice(stat_all, k, replace=False).mean()
                             for _ in range(n_perm)])
            p_two = float((np.abs(null - null.mean())
                           >= abs(obs - null.mean())).mean())
            conf_rows.append({"day": day, "module": name, "n_genes": k,
                              "adj_mean_interaction": obs,
                              "z": float((obs - null.mean())
                                         / (null.std() or np.nan)),
                              "p_perm_adjusted": p_two})

    cdf = pd.DataFrame(conf_rows)
    cdf["q_BH"] = _bh(cdf["p_perm_adjusted"].values)
    cdf = cdf.sort_values(["day", "p_perm_adjusted"])
    cdf.to_csv(OUT / "h3_module_tests_RECOVERY_ADJUSTED.csv", index=False)
    log("\nModule tests AFTER adjusting for recovery and abundance:")
    log("\n" + cdf.to_string(index=False))

    log("\n" + "-" * 74)
    fa = cdf[(cdf.module == "FOCAL")]
    sig_adj = fa[fa.q_BH < cfg.ALPHA]
    log("VERDICT ON THE PRIMARY ENDPOINT")
    if len(sig_adj) == len(fa) and len(fa):
        log(f"  FOCAL survives adjustment at {len(sig_adj)}/{len(fa)} "
            "timepoints. The interaction is NOT explained by read recovery.")
        log("  H3 stands.")
    elif len(sig_adj):
        log(f"  FOCAL survives at only {len(sig_adj)}/{len(fa)} timepoints "
            "after adjustment - partially confounded. Report both.", "WARN")
    else:
        log("  FOCAL does NOT survive adjustment. The unadjusted result is "
            "explained by relative read recovery, exactly as H1's 3'UTR "
            "signal was explained by sequence length. H3 IS NOT SUPPORTED.",
            "ERROR")

    # ------------------------------------------------------------------
    # PER-GENE ATTRIBUTION AND LEAVE-ONE-OUT
    # ------------------------------------------------------------------
    # Two things this is guarding against, both learned the hard way in H1:
    #
    #  1. In H1 the single largest contributor to significance was TGFB1 -
    #     the gene with the WORST sequence recovery. A module-level result
    #     driven by its most poorly-measured member is not a module result.
    #     The recovery adjustment above is linear with a near-zero slope, so
    #     it barely corrects genes in the extreme tail, where log2FC
    #     attenuation is nonlinear. Linear adjustment is therefore NOT
    #     sufficient reassurance for genes at the 0.6th percentile.
    #
    #  2. In H1 the unadjusted result crossed 0.05 in one build and not
    #     another, differing by a single control gene. Leave-one-out shows
    #     directly whether one member carries the result.
    log("\n" + "=" * 74)
    log("PER-GENE ATTRIBUTION - which focal genes carry the result?")
    log("=" * 74)

    fg_rows = []
    for g in sorted(focal_set):
        if g not in res[2].index:
            continue
        fg_rows.append({
            "gene": g,
            "interaction_d2": round(float(res[2].loc[g, "interaction"]), 3),
            "interaction_d5": round(float(res[5].loc[g, "interaction"]), 3),
            "lfc_acomys_d5": round(float(res[5].loc[g, "lfc_acomys"]), 3),
            "lfc_mus_d5": round(float(res[5].loc[g, "lfc_mus"]), 3),
            "tpm_acomys": round(float(aco.loc[g].mean()), 1),
            "tpm_mus": round(float(mus.loc[g].mean()), 1),
            "recovery_pctile": round(float(sstats.percentileofscore(
                recovery.values, recovery.loc[g])), 1),
        })
    fdf = pd.DataFrame(fg_rows)
    fdf.to_csv(OUT / "h3_focal_gene_detail.csv", index=False)
    log("\n" + fdf.to_string(index=False))
    log("\n  recovery_pctile < 5 means the gene is in the bottom 5% for "
        "Acomys recovery relative to Mus - its log2FC is the most likely to "
        "be attenuated by the pseudocount floor.")

    # Leave-one-out is applied to the focal module AND to any panel block that
    # reached significance. A newly-significant block deserves the same
    # scrutiny as the primary endpoint - more so where its membership depends
    # on an ortholog assignment recovered by homology. HAVCR1 in the tubular
    # block cleared its reciprocal-best-hit check by only 1.51x (against
    # TIMD7-PS, with TIMD2 and TIMD6 also hitting), the weakest of the three
    # recovered aliases. If the block's significance rests on that gene, the
    # finding rests on an orthology call rather than on biology.
    loo_targets = {"FOCAL": sorted(focal_set)}
    for name in mods:
        if not name.startswith("PANEL_") or name in ("PANEL_combined15",
                                                     "PANEL_full20"):
            continue
        sig = tab[(tab.module == name) & (tab.q_BH < cfg.ALPHA)]
        if len(sig):
            loo_targets[name] = mods[name]

    for tgt, genes in loo_targets.items():
        log(f"\nLEAVE-ONE-OUT on {tgt}:")
        for day in (2, 5):
            inter = res[day]["interaction"]
            present = [g for g in genes if g in inter.index]
            stat_all = inter.values
            log(f"  day {day}:")
            for drop in present:
                keep = [g for g in present if g != drop]
                if len(keep) < 2:
                    continue
                obs = float(inter.loc[keep].mean())
                k = len(keep)
                null = np.array([rng.choice(stat_all, k, replace=False).mean()
                                 for _ in range(4000)])
                p = float((np.abs(null - null.mean())
                           >= abs(obs - null.mean())).mean())
                flag = "  <-- result depends on this gene" if p > 0.05 else ""
                log(f"    without {drop:8s} n={k}  mean={obs:+.3f}  "
                    f"p={p:.4f}{flag}")

        # Per-gene detail, including recovery percentile, for any block that
        # fired. Genes recovered from LOC placeholders are marked.
        if tgt == "FOCAL":
            continue
        alias = acomys_aliases()
        log(f"\n  per-gene detail for {tgt}:")
        log(f"    {'gene':10s} {'int_d2':>8s} {'int_d5':>8s} "
            f"{'tpm_aco':>9s} {'tpm_mus':>9s} {'recov_%ile':>10s}  origin")
        for g in sorted(genes):
            if g not in res[2].index:
                log(f"    {g:10s} {'absent from shared space':>50s}")
                continue
            origin = (f"recovered from {alias[g]}" if g in alias
                      else "native symbol")
            log(f"    {g:10s} {res[2].loc[g,'interaction']:8.3f} "
                f"{res[5].loc[g,'interaction']:8.3f} "
                f"{float(aco.loc[g].mean()):9.1f} "
                f"{float(mus.loc[g].mean()):9.1f} "
                f"{float(sstats.percentileofscore(recovery.values, recovery.loc[g])):10.1f}"
                f"  {origin}")

    # ------------------------------------------------------------------
    # SELF-CONTAINED TEST - does the module beat SAMPLING noise?
    # ------------------------------------------------------------------
    # The competitive permutation above asks whether the module is shifted
    # relative to other genes. It treats each gene's log2FC as a fixed
    # quantity and never touches replicate variability, so at n=3 it can
    # report a confident p-value for a difference that three animals per
    # group cannot actually resolve. The pre-registration named ROAST as well
    # as CAMERA for exactly this reason.
    #
    # Here the SPECIES LABELS are permuted and the module statistic recomputed.
    #
    # The permutation is STRATIFIED: labels are shuffled within the sham
    # samples and within the injured samples separately. An unstratified
    # shuffle across all 12 would allow a permuted "species" made of 5 shams
    # and 1 injured animal - a configuration the real design can never
    # produce - which distorts the null. Stratifying gives
    # C(6,3) x C(6,3) = 400 label assignments per timepoint, so the finest
    # attainable p is 1/400 = 0.0025.
    from itertools import combinations

    log("\n" + "=" * 74)
    log("SELF-CONTAINED TEST - species labels permuted, stratified by "
        "condition")
    log("=" * 74)
    both = pd.concat([aco, mus], axis=1)
    sc_rows = []
    for day in (2, 5):
        sham_all = [g for g, sp, c, d in samples() if c == "sham"]
        inj_all = [g for g, sp, c, d in samples() if c == "UUO" and d == day]
        true_sham = tuple(g for g in sham_all
                          if meta[g][0] == "Acomys cahirinus")
        true_inj = tuple(g for g in inj_all
                         if meta[g][0] == "Acomys cahirinus")

        def module_stat(a_sham, a_inj):
            m_sham = [c for c in sham_all if c not in a_sham]
            m_inj = [c for c in inj_all if c not in a_inj]
            a = (np.log2(both[list(a_inj)].mean(axis=1) + 1)
                 - np.log2(both[list(a_sham)].mean(axis=1) + 1))
            m = (np.log2(both[m_inj].mean(axis=1) + 1)
                 - np.log2(both[m_sham].mean(axis=1) + 1))
            return a - m

        obs_int = module_stat(true_sham, true_inj)
        # Use the ACTUAL number of Acomys samples in each stratum, not
        # len//2. Under --exclude the groups are no longer balanced, and
        # splitting in half would generate assignments of the wrong shape -
        # the true labelling would not even be among them, so the observed
        # statistic would be compared against a null it cannot belong to.
        assignments = [(s, i)
                       for s in combinations(sham_all, len(true_sham))
                       for i in combinations(inj_all, len(true_inj))]
        cache = {key: module_stat(*key) for key in assignments}

        for name, genes in mods.items():
            present = [g for g in genes if g in both.index]
            if len(present) < 4:
                continue
            obs = float(obs_int.loc[present].mean())
            null = np.array([float(v.loc[present].mean())
                             for v in cache.values()])
            p = float((np.abs(null) >= abs(obs)).mean())
            sc_rows.append({"day": day, "module": name,
                            "n_genes": len(present),
                            "observed": round(obs, 3),
                            "n_permutations": len(null),
                            "p_selfcontained": round(p, 4)})

    if sc_rows:
        sdf = pd.DataFrame(sc_rows).sort_values(["day", "p_selfcontained"])
        sdf.to_csv(OUT / "h3_selfcontained.csv", index=False)
        log("\n" + sdf.to_string(index=False))
        log(f"\nNOTE: {len(assignments)} stratified label assignments, so the "
            f"finest attainable p is {1/len(assignments):.4f}. The true "
            "labelling is itself one of these, so p can never be 0.")

    # SECONDARY - ratio trajectories
    log("\n" + "-" * 74)
    log("SECONDARY - ratio trajectories (log2)")
    pairs = [("TGFB1", "TGFB3"), ("CCR2", "CX3CR1")]
    rrows = []
    for num, den in pairs:
        # Index is uppercase throughout; no case fallback needed.
        if num not in aco.index or den not in aco.index:
            log(f"  {num}:{den} - not both present in shared space", "WARN")
            continue
        for label, df in (("Acomys", aco), ("Mus", mus)):
            for cond, day in (("sham", 0), ("UUO", 2), ("UUO", 5)):
                c = cols_for(df, cond, day)
                r = np.log2((df.loc[num, c].mean() + 1)
                            / (df.loc[den, c].mean() + 1))
                rrows.append({"ratio": f"{num}:{den}", "species": label,
                              "condition": cond, "day": day,
                              "log2_ratio": round(float(r), 3)})
    if rrows:
        rdf = pd.DataFrame(rrows)
        rdf.to_csv(OUT / "h3_ratios.csv", index=False)
        log("\n" + rdf.pivot_table(index=["ratio", "day"], columns="species",
                                   values="log2_ratio").to_string())

    # EXPLORATORY - never a claim
    for day, df in res.items():
        df.sort_values("interaction").to_csv(
            OUT / f"h3_gene_level_EXPLORATORY_day{day}.csv")
    log("\nGene-level tables written with EXPLORATORY in the filename. "
        "n=3 per group; no gene-level claim is made from them.")

    log(f"\nwrote {OUT / 'h3_module_tests.csv'}")
    return True


def _bh(p):
    import numpy as np
    p = np.asarray(p, float)
    n = len(p)
    o = np.argsort(p)
    q = np.empty(n)
    prev = 1.0
    for rank, i in enumerate(o[::-1]):
        prev = min(prev, p[i] * n / (n - rank))
        q[i] = prev
    return q


def main():
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--check", action="store_true")
    ap.add_argument("--runs", action="store_true",
                    help="resolve GSM -> SRR accessions via NCBI")
    ap.add_argument("--force", action="store_true",
                    help="re-query NCBI even if cached")
    ap.add_argument("--fetch-ref", action="store_true",
                    help="download both RefSeq transcriptomes")
    ap.add_argument("--index", action="store_true",
                    help="build kallisto indices")
    ap.add_argument("--quant", action="store_true",
                    help="fetch, quantify and delete each run in turn")
    ap.add_argument("--keep-fastq", action="store_true",
                    help="do not delete FASTQ after quantification "
                         "(needs ~150 GB free)")
    ap.add_argument("--analyse", action="store_true",
                    help="H3 module-level tests")
    ap.add_argument("--symbols", action="store_true",
                    help="diagnose the cross-species symbol namespaces")
    ap.add_argument("--find-gene", metavar="SYMBOL",
                    help="locate a gene in the Acomys transcriptome by "
                         "homology (for LOC-named genes such as CCR2)")
    ap.add_argument("--exclude", metavar="GSM[,GSM...]",
                    help="drop samples, for sensitivity analyses "
                         "(e.g. --exclude GSM5171803)")
    a = ap.parse_args()

    if a.exclude:
        EXCLUDED_SAMPLES.update(x.strip() for x in a.exclude.split(",")
                                if x.strip())
        log("=" * 74, "WARN")
        log(f"SENSITIVITY RUN - excluding {sorted(EXCLUDED_SAMPLES)}", "WARN")
        log("Results below are NOT the primary analysis. Do not mix them "
            "with full-cohort numbers in the same table.", "WARN")
        left = {}
        for _, sp, cond, day in samples():
            left[(sp, cond, day)] = left.get((sp, cond, day), 0) + 1
        for k in sorted(left):
            log(f"  remaining {k[0]:17s} {k[1]:5s} d{k[2]}  n={left[k]}", "WARN")
        if min(left.values()) < 2:
            log("  A group has fewer than 2 replicates - the interaction is "
                "not estimable. Aborting.", "ERROR")
            sys.exit(1)
        log("=" * 74, "WARN")
    if a.check:
        sys.exit(0 if check() else 1)
    elif a.runs:
        resolve_runs(force=a.force)
    elif a.fetch_ref:
        sys.exit(0 if fetch_reference() else 1)
    elif a.index:
        sys.exit(0 if build_indices() else 1)
    elif a.quant:
        sys.exit(0 if quant_all(keep_fastq=a.keep_fastq) else 1)
    elif a.analyse:
        sys.exit(0 if analyse() else 1)
    elif a.symbols:
        sys.exit(0 if diagnose_symbols() else 1)
    elif a.find_gene:
        sys.exit(0 if find_gene(a.find_gene) else 1)
    else:
        ap.print_help()


if __name__ == "__main__":
    main()
