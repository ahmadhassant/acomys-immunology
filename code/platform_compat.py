#!/usr/bin/env python3
"""
================================================================================
PLATFORM COMPATIBILITY - native Windows and Linux from one codebase
================================================================================

The pipeline runs natively on Windows. No WSL, no conda, no Linux VM.

TOOL SUBSTITUTIONS (all official, all Windows-native)

    purpose              Linux            Windows-native        note
    -------------------  ---------------  --------------------  ----------------
    genome download      datasets         datasets.exe          same tool
    SRA download         sra-tools        SRA Toolkit win64     same tool
    homology / orthology BLAST+           BLAST+ win64          same tool
    alignment            MAFFT / PRANK    MAFFT for Windows     same tool
    selection tests      HyPhy aBSREL     PAML codeml           branch-site
                                                                model A test
    quantification       Salmon           kallisto (win64)      pseudoalignment
    differential expr.   R limma/edgeR    PyDESeq2              pure Python
    UTR evidence         HISAT2+StringTie congener annotation   see below

The last row is not a downgrade. Rather than aligning 105M pooled reads to
infer UTR boundaries from coverage, we take them from *Acomys dimidiatus*,
whose Ensembl annotation is manually curated and which sits ~2 Ma from the
focal species - closer than *A. russatus* at ~5 Ma. Curated congener
boundaries beat coverage-inferred ones from pooled tissue, and it removes a
260 GB download and the entire alignment step.

THE WINDOWS PATH PROBLEM
------------------------
BLAST parses -in / -query / -db as SPACE-SEPARATED FILE LISTS, so a path like
    C:\\my research\\shared papers\\...
is read as several files and fails. On Linux we stage symlinks; on Windows
symlinks need admin rights or Developer Mode, so we use the 8.3 short path
instead (C:\\MYRESE~1\\SHARED~1\\...), which is space-free by construction and
costs nothing.

Usage:
    import platform_compat as pc
    pc.report()                     # what is installed, what is missing
    subprocess.run([pc.tool("blastn"), "-query", pc.safe_path(q), ...])
================================================================================
"""

from __future__ import annotations

import ctypes
import gzip
import os
import platform
import shutil
import subprocess
from pathlib import Path

IS_WINDOWS = platform.system() == "Windows"
IS_LINUX = platform.system() == "Linux"


# ==============================================================================
# PATHS
# ==============================================================================

def short_path(p: Path | str) -> str:
    """Windows 8.3 short path - space-free, no admin rights needed."""
    p = str(p)
    if not IS_WINDOWS:
        return p
    try:
        buf = ctypes.create_unicode_buffer(len(p) + 260)
        n = ctypes.windll.kernel32.GetShortPathNameW(p, buf, len(buf))
        if n and n < len(buf) and " " not in buf.value:
            return buf.value
    except Exception:
        pass
    return p


def safe_path(p: Path | str, scratch: Path | None = None) -> str:
    """A path BLAST can parse: space-free, and still pointing at `p`.

    Handles three cases that all arise in this pipeline:
      - an existing FILE      (a query FASTA, a genome)
      - an existing DIRECTORY (an output folder)
      - a path that does NOT exist yet (a BLAST db PREFIX like .../GCA_x_db,
        which is never a real file - only <prefix>.nin etc. are)

    Windows -> 8.3 short path. If the target does not exist, shorten its
               PARENT and re-attach the name, which is what makes db prefixes
               work.
    Linux   -> symlink staged under scratch (files only).
    """
    p = Path(p)
    s = str(p)
    if " " not in s:
        return s

    if IS_WINDOWS:
        # 1. direct short path (works when p exists)
        sp = short_path(p)
        if " " not in sp:
            return sp
        # 2. shorten the parent and re-attach the name. Covers non-existent
        #    paths such as BLAST database prefixes.
        parent = short_path(p.parent)
        if " " not in parent:
            return str(Path(parent) / p.name)
        # 3. 8.3 generation is disabled on this volume. Copy, but only a
        #    regular file - never a directory.
        if p.is_file():
            stage = Path(scratch or os.environ.get("ACOMYS_SCRATCH")
                         or Path.home() / "acomys_scratch") / "blast_stage"
            stage.mkdir(parents=True, exist_ok=True)
            dest = stage / p.name.replace(" ", "_")
            if not dest.exists() or dest.stat().st_size != p.stat().st_size:
                shutil.copy2(p, dest)
            return str(dest)
        # 4. nothing worked - the caller should relocate the project.
        return s

    import hashlib
    if not p.is_file():
        return s
    stage = Path(scratch or os.environ.get("ACOMYS_SCRATCH") or "/tmp") / "blast_stage"
    stage.mkdir(parents=True, exist_ok=True)
    link = stage / f"{hashlib.md5(s.encode()).hexdigest()[:10]}_{p.name.replace(' ', '_')}"
    try:
        if not link.exists():
            link.symlink_to(p)
        return str(link)
    except OSError:
        return f'"{s}"'


def eight_dot_three_available(path: Path | str) -> bool:
    """Does this volume generate 8.3 short names? Decides whether a project
    path containing spaces is workable at all."""
    if not IS_WINDOWS:
        return False
    p = Path(path)
    probe = p if p.exists() else p.parent
    sp = short_path(probe)
    return " " not in sp and sp != str(probe)


def space_free_dir(preferred: Path) -> Path:
    """A writable directory with no spaces, for tool outputs."""
    preferred = Path(preferred)
    if " " not in str(preferred):
        preferred.mkdir(parents=True, exist_ok=True)
        return preferred
    root = Path(os.environ.get("ACOMYS_SCRATCH")
                or (Path.home() / "acomys_scratch")) / preferred.name
    root.mkdir(parents=True, exist_ok=True)
    return root


# ==============================================================================
# TOOL DISCOVERY
# ==============================================================================
# Windows installers do not add themselves to PATH, so search the usual places.

# Patterns use a leading '*:' which is expanded across all fixed drives, because
# the project has lived on both C: and D:. Set ACOMYS_TOOL_DIRS (os.pathsep
# separated) to prepend extra search directories.
_WIN_HINTS = {
    "blastn":      [r"*:\Program Files\NCBI\blast*\bin", r"*:\NCBI\blast*\bin"],
    "makeblastdb": [r"*:\Program Files\NCBI\blast*\bin", r"*:\NCBI\blast*\bin"],
    "blastdbcmd":  [r"*:\Program Files\NCBI\blast*\bin", r"*:\NCBI\blast*\bin"],
    "datasets":    [r"*:\Tools\ncbi-datasets", r"*:\Program Files\NCBI\datasets"],
    "prefetch":    [r"*:\Tools\sratoolkit*\bin", r"*:\Program Files\sratoolkit*\bin"],
    "fasterq-dump": [r"*:\Tools\sratoolkit*\bin", r"*:\Program Files\sratoolkit*\bin"],
    # The Windows archive extracts to a versioned folder
    # (kallisto_windows-v0.51.1\), sometimes with a nested kallisto\ inside,
    # so glob rather than matching an exact directory name.
    "kallisto":    [r"*:\Tools\kallisto*", r"*:\Tools\kallisto*\kallisto",
                    r"*:\Tools\*\kallisto", r"*:\Program Files\kallisto*"],
    # PAML win release extracts as paml-<ver>-win-x86_64\bin\codeml.exe
    "codeml":      [r"*:\Tools\paml*\bin", r"*:\paml*\bin", r"*:\Tools\paml*"],
    # MAFFT all-in-one extracts as mafft-win\mafft.bat (NOT an .exe).
    "mafft":       [r"*:\Tools\mafft*", r"*:\Tools\mafft*\mafft-win",
                    r"*:\Program Files\mafft*", r"*:\mafft-win"],
}

# Windows executable suffixes, in preference order. '.bat' matters: the MAFFT
# all-in-one package ships mafft.bat wrapping an MSYS tree, with no mafft.exe.
_WIN_EXTS = (".exe", ".bat", ".cmd")


def _expand_drives(pattern: str) -> list[str]:
    """Expand a leading '*:' into every fixed drive present on the machine."""
    if not pattern.startswith("*:"):
        return [pattern]
    out = []
    for letter in "CDEFGH":
        root = f"{letter}:\\"
        if Path(root).exists():
            out.append(letter + ":" + pattern[2:])
    return out

_cache: dict[str, str | None] = {}


def tool(name: str) -> str | None:
    """Absolute path to an external tool, or None. Cached.

    On Windows this resolves .exe, .bat and .cmd. A .bat result is NOT
    directly executable by subprocess - always launch via run_tool().
    """
    if name in _cache:
        return _cache[name]

    found = shutil.which(name)

    if not found and IS_WINDOWS:
        for ext in _WIN_EXTS:
            found = shutil.which(name + ext)
            if found:
                break

    if not found and IS_WINDOWS:
        import glob
        # User-supplied directories win over the built-in hints.
        extra = [d for d in os.environ.get("ACOMYS_TOOL_DIRS", "").split(os.pathsep) if d]
        patterns = extra + _WIN_HINTS.get(name, [])
        for pattern in patterns:
            for expanded in _expand_drives(pattern):
                for d in glob.glob(expanded):
                    for ext in _WIN_EXTS:
                        cand = Path(d) / (name + ext)
                        if cand.exists():
                            found = str(cand)
                            break
                    if found:
                        break
                if found:
                    break
            if found:
                break

    _cache[name] = found
    return found


def run_tool(exe: str, args: list[str], **kw):
    """subprocess.run for an external tool, safe for Windows .bat wrappers.

    CreateProcess cannot execute a .bat directly - Python raises
    'WinError 193: %1 is not a valid Win32 application'. Batch files must be
    launched through the command interpreter. This is exactly how the MAFFT
    all-in-one Windows package is distributed, so calling mafft without this
    wrapper fails at the first alignment.
    """
    import subprocess
    cmd = [exe, *args]
    if IS_WINDOWS and exe.lower().endswith((".bat", ".cmd")):
        # /c = run then terminate. The path is passed as its own argv entry so
        # spaces in it are handled by cmd's own quoting, not ours.
        cmd = ["cmd.exe", "/c", exe, *args]
    kw.setdefault("check", True)
    return subprocess.run(cmd, **kw)


def require(name: str) -> str:
    p = tool(name)
    if p is None:
        raise FileNotFoundError(
            f"'{name}' not found. On Windows run setup_windows.ps1; "
            f"on Linux create the conda envs. See platform_compat.report().")
    return p


# ==============================================================================
# PORTABLE UTILITIES
# ==============================================================================

def gzip_file(path: Path, remove_original: bool = True) -> Path:
    """Compress in pure Python - no external gzip, which Windows lacks."""
    path = Path(path)
    dest = path.with_suffix(path.suffix + ".gz")
    with open(path, "rb") as fin, gzip.open(dest, "wb", compresslevel=6) as fout:
        shutil.copyfileobj(fin, fout, length=1 << 20)
    if remove_original:
        path.unlink()
    return dest


def unzip(archive: Path, dest: Path, members_endswith: str | None = None):
    """Extract with zipfile - no external unzip binary required."""
    import zipfile
    dest = Path(dest)
    dest.mkdir(parents=True, exist_ok=True)
    with zipfile.ZipFile(archive) as zf:
        names = zf.namelist()
        if members_endswith:
            names = [n for n in names if n.endswith(members_endswith)]
        for n in names:
            zf.extract(n, dest)


def run(cmd: list, **kw) -> subprocess.CompletedProcess:
    """subprocess.run with Windows-friendly defaults."""
    kw.setdefault("capture_output", True)
    if IS_WINDOWS:
        # Stop console windows flashing up for every BLAST call.
        kw.setdefault("creationflags", getattr(subprocess, "CREATE_NO_WINDOW", 0))
    return subprocess.run([str(c) for c in cmd], **kw)


# ==============================================================================
# REPORT
# ==============================================================================

REQUIRED = {
    "blastn":       "WP2 orthology",
    "makeblastdb":  "WP2 orthology",
    "blastdbcmd":   "WP2 orthology (optional - sseq route preferred)",
    "datasets":     "WP1 genome download",
}
OPTIONAL = {
    "prefetch":     "WP5 raw reads",
    "fasterq-dump": "WP5 raw reads",
    "kallisto":     "WP5 quantification (replaces Salmon)",
    "codeml":       "WP4 selection tests (replaces HyPhy)",
    "mafft":        "WP4 alignment",
}


def report() -> bool:
    print("=" * 70)
    print(f"PLATFORM: {platform.system()} {platform.release()}")
    print(f"python  : {platform.python_version()}")
    print("=" * 70)
    ok = True
    print("\nRequired:")
    for t, why in REQUIRED.items():
        p = tool(t)
        print(f"  [{'OK     ' if p else 'MISSING'}] {t:14s} {why}")
        if not p and "optional" not in why:
            ok = False
    print("\nOptional (later work packages):")
    for t, why in OPTIONAL.items():
        p = tool(t)
        print(f"  [{'OK     ' if p else 'missing'}] {t:14s} {why}")

    scratch = os.environ.get("ACOMYS_SCRATCH")
    print(f"\nACOMYS_SCRATCH: {scratch or 'NOT SET'}")

    root = Path(__file__).resolve().parent.parent
    if " " in str(root):
        print(f"\nProject path contains a space:\n  {root}")
        if eight_dot_three_available(root):
            print(f"  8.3 short path : {short_path(root)}")
            print("  Handled automatically - BLAST will use the short form.")
        else:
            ok = False
            print("  8.3 SHORT NAMES ARE DISABLED ON THIS VOLUME.")
            print("  BLAST parses -in/-query as space-separated file LISTS, so")
            print("  it cannot open a path containing a space, and there is no")
            print("  workaround left on this drive.")
            print()
            print("  FIX - move the project somewhere without spaces, e.g.:")
            print(f"    move \"{root}\" D:\\acomys")
            print("  (or enable 8.3 as admin: fsutil 8dot3name set D: 0 -")
            print("   note this only affects folders created afterwards)")

    print("\n" + "=" * 70)
    print(f"READY: {ok}")
    if not ok:
        print("Run  powershell -ExecutionPolicy Bypass -File setup_windows.ps1")
    print("=" * 70)
    return ok


if __name__ == "__main__":
    report()
