#!/usr/bin/env python3
"""
================================================================================
RUN ALL - single entry point for the whole pipeline
Acomys cahirinus regeneration project
================================================================================

Follows the numbered-step pattern from FISH/new update/analysis/ - the most
reproducible thing in the papers folder.

    python run_all.py --check          # environment + data readiness, no work
    python run_all.py --dry-run        # show what each step would do
    python run_all.py --from 2         # resume from step 2
    python run_all.py --only 3         # one step
    python run_all.py                  # everything

Each step is idempotent: re-running skips work already done (cached downloads,
existing outputs). Nothing is ever edited in place.
================================================================================
"""

from __future__ import annotations

import argparse
import os
import subprocess
import sys
import time
from pathlib import Path

CODE = Path(__file__).resolve().parent / "code"
sys.path.insert(0, str(CODE))
import config as cfg

log = cfg.log

STEPS = [
    (1, "01_fetch_data.py", ["--all"],
     "WP1  Acquire GEO tables, genomes, transcriptome, orthologs"),
    (2, "02_build_ortholog_scaffold.py", ["--all"],
     "WP2  Build 1:1 ortholog scaffold + region slicing  [LINCHPIN]"),
    (3, "03_kmer_phylo_controlled.py", ["--arm", "all", "--k-sweep"],
     "WP3  Phylogeny-controlled alignment-free analysis  [PRIMARY]"),
    # (4, "04_selection_analysis.py", [], "WP4  HyPhy aBSREL / RELAX"),
    # (5, "05_expression_reanalysis.py", [], "WP5  Species-native re-quantification"),
    # (6, "06_integration.py", [], "WP6  Sequence x expression convergence"),
]


SUPPORTED_PY = (3, 10, 3, 12)   # min_major, min_minor, max_major, max_minor

CONDA_NAME = {"sklearn": "scikit-learn", "Bio": "biopython",
              "matplotlib": "matplotlib-base"}


def preflight() -> bool:
    log("=" * 74)
    log("PREFLIGHT")
    log("=" * 74)
    ok = True

    # ---- is `python` actually the ACTIVE env's python? -------------------
    # The commonest silent failure: `conda env create` fails at solve time,
    # leaving an env with no interpreter. `conda activate` still succeeds, and
    # `python` falls through to base - so you are running a different env than
    # you think you are.
    conda_prefix = os.environ.get("CONDA_PREFIX")
    conda_name = os.environ.get("CONDA_DEFAULT_ENV", "?")
    log(f"python {sys.version.split()[0]}   ({sys.executable})")
    if conda_prefix:
        log(f"active conda env: {conda_name}  ({conda_prefix})")
        # Compare RESOLVED PATHS, not string prefixes: "/envs/acomys-bio"
        # startswith "/envs/acomys" is True, which would silently mask the
        # exact shadowing case this check exists to catch.
        if Path(sys.prefix).resolve() != Path(conda_prefix).resolve():
            ok = False
            log("", "ERROR")
            log("  MISMATCH: the running python is NOT from the active env.",
                "ERROR")
            log(f"    active env : {conda_prefix}", "ERROR")
            log(f"    python from: {sys.prefix}", "ERROR")

            envs_dir = str(Path(conda_prefix).parent)
            if Path(envs_dir) in Path(sys.prefix).resolve().parents:
                # Running another conda env's python -> PATH shadowing.
                other = Path(sys.prefix).name
                log("", "ERROR")
                log(f"  CAUSE: PATH shadowing. Env '{other}' was PREPENDED to "
                    f"PATH, and it ships its own python, so it wins over the "
                    f"activated env. Every package check below describes "
                    f"'{other}', not '{conda_name}'.", "ERROR")
                log("", "ERROR")
                log("  Never PREPEND a conda env's bin/ to PATH - conda envs "
                    "contain interpreters, not just tools.", "ERROR")
                log("  FIX - expose only the binaries you need, never the "
                    "interpreter:", "ERROR")
                log(f"    for t in datasets blastn blastp makeblastdb diamond \\",
                    "ERROR")
                log(f"             orthofinder hyphy mafft prank salmon \\", "ERROR")
                log(f"             fastqc multiqc prefetch fasterq-dump; do", "ERROR")
                log(f"      ln -sf {envs_dir}/{other}/bin/$t \\", "ERROR")
                log(f"             {conda_prefix}/bin/$t 2>/dev/null", "ERROR")
                log(f"    done", "ERROR")
                log("  Then remove the bad export from your shell and "
                    "~/.bashrc, and re-activate.", "ERROR")
                log("  (If you prefer PATH, APPEND it: "
                    'export PATH="$PATH:'f'{envs_dir}/{other}/bin")', "ERROR")
            else:
                log("  CAUSE: env creation likely FAILED and left no "
                    "interpreter, so `python` fell through to base.", "ERROR")
                log("  Diagnose:", "ERROR")
                log("    conda env list", "ERROR")
                log("    ls $CONDA_PREFIX/bin/python*", "ERROR")
                log("  Rebuild step by step so failures are visible:", "ERROR")
                log("    conda deactivate && conda env remove -n acomys -y", "ERROR")
                log("    conda create -n acomys python=3.11 -y", "ERROR")
                log("    conda activate acomys", "ERROR")
                log("    conda install -c conda-forge numpy pandas scipy "
                    "scikit-learn statsmodels matplotlib-base seaborn openpyxl "
                    "requests tqdm biopython -y", "ERROR")
                log("    pip install shap dendropy", "ERROR")

    # ---- Python version -------------------------------------------------
    v = sys.version_info
    lo = (SUPPORTED_PY[0], SUPPORTED_PY[1])
    hi = (SUPPORTED_PY[2], SUPPORTED_PY[3])
    if not (lo <= (v.major, v.minor) <= hi):
        ok = False
        log(f"  UNSUPPORTED Python. Need {lo[0]}.{lo[1]}-{hi[0]}.{hi[1]}; "
            f"environment.yml pins 3.11.", "ERROR")
        if (v.major, v.minor) > hi:
            log("  Python is NEWER than the scientific stack supports. conda "
                "will silently drop packages that have no build for it - which "
                "is what a half-populated env below means.", "ERROR")
        log("  This env was almost certainly not built from environment.yml.",
            "ERROR")
        log("  Rebuild:", "ERROR")
        log("    conda deactivate", "ERROR")
        log("    conda env remove -n acomys -y", "ERROR")
        log("    conda env create -f environment.yml", "ERROR")
        log("    conda activate acomys", "ERROR")

    # ---- packages -------------------------------------------------------
    missing = []
    for mod in ["numpy", "pandas", "scipy", "sklearn", "Bio",
                "matplotlib", "seaborn"]:
        try:
            __import__(mod)
            log(f"  [OK     ] {mod}")
        except ImportError:
            log(f"  [MISSING] {mod}")
            missing.append(CONDA_NAME.get(mod, mod))
            ok = False
    for mod in ["shap", "statsmodels"]:
        try:
            __import__(mod)
            log(f"  [OK     ] {mod}")
        except ImportError:
            log(f"  [optional] {mod}")

    if missing:
        log("")
        log(f"  Install the {len(missing)} missing package(s):", "ERROR")
        log(f"    conda install -c conda-forge {' '.join(missing)}", "ERROR")
        log("  If several are missing at once, rebuilding from "
            "environment.yml is safer than patching one by one.", "ERROR")

    # ---- external CLI tools ---------------------------------------------
    from shutil import which
    log("")
    log("External tools (environment-bio.yml):")
    for tool, needed_by in [("datasets", "WP1 genomes"),
                            ("blastn", "WP2 orthology"),
                            ("makeblastdb", "WP2 orthology"),
                            ("salmon", "WP5 quantification"),
                            ("hyphy", "WP4 selection")]:
        present = which(tool) is not None
        log(f"  [{'OK     ' if present else 'MISSING'}] {tool:12s} {needed_by}")
    if not which("datasets") or not which("blastn"):
        log("  Symlink the binaries into this env. Do NOT prepend "
            "acomys-bio/bin to PATH - it ships its own python and will "
            "shadow the interpreter you just activated:", "WARN")
        log("      for t in datasets blastn makeblastdb diamond orthofinder \\", "WARN")
        log("               hyphy mafft prank salmon fastqc multiqc; do", "WARN")
        log("        ln -sf $(conda info --base)/envs/acomys-bio/bin/$t \\", "WARN")
        log("               $CONDA_PREFIX/bin/$t 2>/dev/null; done",
            "WARN")

    log("")
    log(f"NCBI_API_KEY: {'set' if os.environ.get('NCBI_API_KEY') else 'NOT SET (3x slower)'}")

    log("")
    cfg.summary()
    if len(cfg.USABLE_SPECIES) < cfg.MIN_SPECIES_REQUIRED:
        ok = False

    log("")
    log("Pre-specified decision rule for H1:")
    for line in cfg._wrap(cfg.H1_DECISION_RULE, 68):
        log(f"  {line}")
    log("=" * 74)
    log(f"PREFLIGHT: {'PASS' if ok else 'FAIL'}")
    log("=" * 74)
    return ok


def run_step(num: int, script: str, args: list[str], desc: str,
             dry_run: bool) -> bool:
    log("")
    log("#" * 74)
    log(f"# STEP {num}: {desc}")
    log(f"# {script} {' '.join(args)}")
    log("#" * 74)
    if dry_run:
        log("DRY RUN - not executed")
        return True

    path = CODE / script
    if not path.exists():
        log(f"NOT YET WRITTEN: {script} - skipping", "WARN")
        return True

    t0 = time.time()
    try:
        subprocess.run([sys.executable, str(path), *args],
                       check=True, cwd=str(CODE))
        log(f"STEP {num} completed in {time.time() - t0:.1f}s")
        return True
    except subprocess.CalledProcessError as e:
        log(f"STEP {num} FAILED (exit {e.returncode})", "ERROR")
        return False


def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--check", action="store_true", help="preflight only")
    ap.add_argument("--dry-run", action="store_true")
    ap.add_argument("--from", dest="start", type=int, default=1)
    ap.add_argument("--only", type=int)
    ap.add_argument("--continue-on-error", action="store_true")
    args = ap.parse_args()

    cfg.ensure_dirs()
    ok = preflight()
    if args.check:
        sys.exit(0 if ok else 1)
    if not ok and not args.continue_on_error:
        log("Preflight failed. Fix the above, or pass --continue-on-error.",
            "ERROR")
        sys.exit(1)

    todo = [s for s in STEPS
            if (args.only is None and s[0] >= args.start) or s[0] == args.only]
    if not todo:
        log("Nothing to run.", "WARN")
        return

    t0 = time.time()
    for num, script, sargs, desc in todo:
        if not run_step(num, script, sargs, desc, args.dry_run):
            if not args.continue_on_error:
                sys.exit(1)

    log("")
    log("=" * 74)
    log(f"PIPELINE FINISHED in {time.time() - t0:.1f}s")
    log(f"Results : {cfg.RESULTS}")
    log(f"Manifest: {cfg.MANIFEST}")
    log("=" * 74)


if __name__ == "__main__":
    main()
