#!/usr/bin/env bash
# =============================================================================
# MIGRATE THE ACOMYS PROJECT TO A LINUX SERVER
# =============================================================================
#
# Run this ON THE SERVER, after copying the project folder across.
#
#   1. On the laptop, copy ONLY the small, hard-won files (see below).
#   2. On the server:  bash migrate_to_server.sh
#
# WHAT TO COPY (a few hundred MB):
#   code/                              the pipeline
#   *.md  *.docx  *.yml                docs, plan, environments
#   .env                               your NCBI key (never committed)
#   data/reference/orthologs/          THE ORTHOLOG SCAFFOLD - seven rounds of
#                                      QC went into this; do not regenerate it
#                                      casually
#   data/provenance_manifest.json      hashes of everything downloaded so far
#
# WHAT NOT TO COPY (~150 GB, faster to re-download than to transfer):
#   data/raw/sra/                      67 GB of .sra
#   data/reference/genomes/            10 GB of assemblies
#   fasterq.tmp.*                      delete these on the laptop, don't copy
#
# Example from the laptop (WSL):
#   rsync -avz --progress \
#     --exclude 'data/raw' --exclude 'data/reference/genomes' \
#     --exclude 'fasterq.tmp.*' --exclude '__pycache__' \
#     "/mnt/c/my research/shared papers/Immunology project/" \
#     user@server:/path/to/acomys/
#
# NOTE the trailing slashes, and QUOTE the source: the laptop path contains a
# space, which is what broke BLAST earlier.
# =============================================================================

set -euo pipefail

PROJECT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$PROJECT_DIR"

echo "=============================================================="
echo "Acomys project - server setup"
echo "project : $PROJECT_DIR"
echo "host    : $(hostname)"
echo "=============================================================="

# --- 0. sanity: no spaces in the path ----------------------------------------
if [[ "$PROJECT_DIR" == *" "* ]]; then
    echo "WARNING: project path contains a space:"
    echo "  $PROJECT_DIR"
    echo "BLAST parses -in/-query as space-separated file LISTS, so this breaks"
    echo "silently. The code stages symlinks to work around it, but a path"
    echo "without spaces is simpler. Consider moving the project."
    echo
fi

# --- 1. resources -------------------------------------------------------------
echo "--- resources ---"
CORES=$(nproc)
MEM_GB=$(free -g | awk '/^Mem:/{print $2}')
echo "cores : $CORES"
echo "memory: ${MEM_GB} GB"
echo "disk  : $(df -h . | awk 'NR==2{print $4}') free here"
echo
echo "Guidance for this pipeline:"
echo "  WP2 BLAST      : 4-8 cores,  8 GB"
echo "  WP2b HISAT2    : 8-16 cores, 16 GB, ~200 GB scratch"
echo "  WP4 HyPhy      : 16+ cores  (embarrassingly parallel across genes)"
echo "  WP5 Salmon     : 8-16 cores, 16 GB"
echo "  NO GPU IS USED ANYWHERE. Do not request a GPU node."
echo

# --- 2. scratch ---------------------------------------------------------------
echo "--- scratch space ---"
DEFAULT_SCRATCH="${HOME}/acomys_scratch"
SCRATCH="${ACOMYS_SCRATCH:-$DEFAULT_SCRATCH}"
mkdir -p "$SCRATCH"/{blastdb,blast_stage,sra,sra_tmp,align}
SCRATCH_FREE=$(df -BG "$SCRATCH" | awk 'NR==2{gsub("G","",$4); print $4}')
echo "scratch: $SCRATCH  (${SCRATCH_FREE} GB free)"
if (( SCRATCH_FREE < 250 )); then
    echo "WARNING: under 250 GB. SRA extraction plus alignment needs roughly"
    echo "that much. Point ACOMYS_SCRATCH at a larger local filesystem."
fi
echo

# --- 3. conda environments ----------------------------------------------------
echo "--- conda ---"
if ! command -v conda >/dev/null 2>&1; then
    echo "conda not found. Install miniforge:"
    echo "  wget https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh"
    echo "  bash Miniforge3-Linux-x86_64.sh -b -p \$HOME/miniforge3"
    echo "  \$HOME/miniforge3/bin/conda init bash && exec bash"
    exit 1
fi
echo "conda: $(conda --version)  base: $(conda info --base)"

for env_file in environment.yml environment-bio.yml; do
    env_name=$(awk '/^name:/{print $2; exit}' "$env_file")
    if conda env list | grep -qE "^${env_name}\s"; then
        echo "  [exists] $env_name"
    else
        echo "  [create] $env_name  from $env_file"
        conda env create -f "$env_file"
    fi
done
echo

# --- 4. expose the bio tools without shadowing python -------------------------
# NEVER prepend acomys-bio/bin to PATH: it ships its own python and will
# hijack the analysis interpreter. Symlink the binaries instead.
echo "--- linking bioinformatics tools ---"
CONDA_BASE=$(conda info --base)
ANALYSIS_BIN="$CONDA_BASE/envs/acomys/bin"
BIO_BIN="$CONDA_BASE/envs/acomys-bio/bin"
linked=0
for t in datasets blastn blastp makeblastdb blastdbcmd diamond orthofinder \
         hyphy mafft prank salmon fastqc multiqc prefetch fasterq-dump \
         hisat2 hisat2-build stringtie samtools fastp gffread; do
    if [[ -x "$BIO_BIN/$t" ]]; then
        ln -sf "$BIO_BIN/$t" "$ANALYSIS_BIN/$t"
        linked=$((linked+1))
    fi
done
echo "  linked $linked tools into $ANALYSIS_BIN"
echo

# --- 5. persist environment ---------------------------------------------------
echo "--- shell setup ---"
MARK="# >>> acomys project >>>"
if ! grep -qF "$MARK" "$HOME/.bashrc" 2>/dev/null; then
    {
        echo ""
        echo "$MARK"
        echo "export ACOMYS_SCRATCH=\"$SCRATCH\""
        echo "# NCBI_API_KEY is read from the project's .env - do not put it here"
        echo "# <<< acomys project <<<"
    } >> "$HOME/.bashrc"
    echo "  appended ACOMYS_SCRATCH to ~/.bashrc"
else
    echo "  ~/.bashrc already configured"
fi
echo

# --- 6. verify ----------------------------------------------------------------
echo "--- verification ---"
echo "Run these now:"
echo
echo "  conda activate acomys"
echo "  export ACOMYS_SCRATCH=\"$SCRATCH\""
echo "  python run_all.py --check"
echo
echo "Then re-fetch the bulk data (not transferred):"
echo "  python code/01_fetch_data.py --step genomes        # ~10 GB"
echo "  python code/01_fetch_data.py --step transcriptome  # ~67 GB .sra"
echo
echo "Long jobs: use tmux so a dropped remote-desktop session cannot kill them."
echo "  tmux new -s acomys"
echo "  python code/01_fetch_data.py --step transcriptome"
echo "  # detach: Ctrl-b then d      reattach: tmux attach -t acomys"
echo
echo "=============================================================="
echo "setup complete"
echo "=============================================================="
