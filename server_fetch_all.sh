#!/usr/bin/env bash
# =============================================================================
# ACOMYS PROJECT - FETCH EVERYTHING FROM NCBI, ON THE SERVER, FROM SCRATCH
# =============================================================================
#
# Self-contained. Needs nothing from any other machine - no project code, no
# uploads. Paste this file onto the server and run it.
#
#   chmod +x server_fetch_all.sh
#   tmux new -s fetch                 # ESSENTIAL: survives a dropped session
#   ./server_fetch_all.sh --all
#   # detach with Ctrl-b then d ; reattach with: tmux attach -t fetch
#
# STAGES (each resumable and independently runnable):
#   setup        miniforge + two conda environments
#   genomes      12 murid/gerbil/cricetid assemblies       ~12 GB
#   geo          GSE168876 processed tables + metadata     ~2 MB
#   sra-utr      PRJNA342864, 2 runs -> FASTQ              ~67 GB -> ~190 GB
#   sra-kidney   SRP310563, 18 runs -> FASTQ               ~90 GB -> ~250 GB
#
#   ./server_fetch_all.sh --stage genomes
#   ./server_fetch_all.sh --all --skip sra-kidney
#   ./server_fetch_all.sh --all --keep-sra        # don't delete .sra after dump
#
# DISK: budget ~600 GB. Check with --stage preflight before committing.
# GPU:  none used anywhere in this project. Do not request a GPU node.
#
# Everything downloaded is SHA256-hashed into fetch_manifest.tsv.
# =============================================================================

set -uo pipefail

# --- configuration -----------------------------------------------------------
PROJECT="${ACOMYS_PROJECT:-$HOME/acomys}"
SCRATCH="${ACOMYS_SCRATCH:-$HOME/acomys_scratch}"
NCBI_API_KEY="${NCBI_API_KEY:-c4724f3af36e5d485f3a24ba94d4db539008}"
EMAIL="ahmad.hassanat@gmail.com"
THREADS="${THREADS:-$(nproc)}"
DUMP_THREADS=$(( THREADS > 6 ? 6 : THREADS ))   # more threads = more memory
KEEP_SRA=0

DATA="$PROJECT/data"
RAW="$DATA/raw"
REF="$DATA/reference"
LOGS="$PROJECT/logs"
MANIFEST="$DATA/fetch_manifest.tsv"

# --- the 12 verified assemblies ----------------------------------------------
# Confirmed against NCBI 13-14 Aug 2026. Clade / regenerative status per
# Riddell et al. 2025 PNAS 122:e2420726122 (only deomyines regenerate).
read -r -d '' GENOMES <<'EOF'
GCA_029890205.1	Acomys_cahirinus	Deomyinae	FOCAL_regenerative
GCF_903995435.1	Acomys_russatus	Deomyinae	regenerative_congener
GCA_907164435.1	Acomys_dimidiatus	Deomyinae	regenerative
GCF_030254825.1	Meriones_unguiculatus	Gerbillinae	sister_clade
GCF_907164565.1	Psammomys_obesus	Gerbillinae	sister_clade
GCF_000001635.27	Mus_musculus	Murinae	comparator_has_injury_data
GCF_015227675.2	Rattus_norvegicus	Murinae	comparator
GCF_947179515.1	Apodemus_sylvaticus	Murinae	comparator
GCF_004785775.1	Grammomys_surdaster	Murinae	comparator
GCF_008632895.1	Mastomys_coucha	Murinae	comparator
GCF_049852395.1	Peromyscus_maniculatus	Cricetidae	outgroup
GCF_017639785.1	Mesocricetus_auratus	Cricetidae	outgroup
EOF

GEO_FILES=(
  "GSE168876_Acomys_day_2_vs_acomys_sham.xlsx"
  "GSE168876_Acomys_day_5_vs_acomys_sham.xlsx"
  "GSE168876_Mouse_day_2_vs_mouse_sham.xlsx"
  "GSE168876_Mouse_day_5_vs_mouse_sham.xlsx"
)
GEO_MATRICES=( "GSE168876-GPL24247_series_matrix.txt.gz"
               "GSE168876-GPL29848_series_matrix.txt.gz" )
UTR_RUNS=( SRR4279903 SRR4279904 )        # PRJNA342864, pooled 15 organs
KIDNEY_PROJECT="PRJNA714406"              # SRP310563, 18 kidney samples

# --- helpers ------------------------------------------------------------------
c_hdr(){ printf '\n\033[1;36m=== %s ===\033[0m\n' "$*"; }
ok(){   printf '\033[0;32m  [ok]\033[0m %s\n' "$*"; }
warn(){ printf '\033[0;33m  [warn]\033[0m %s\n' "$*"; }
err(){  printf '\033[0;31m  [ERROR]\033[0m %s\n' "$*"; }
info(){ printf '  %s\n' "$*"; }

record(){   # record <file> <source> <note>
  [[ -f "$1" ]] || return 0
  local h; h=$(sha256sum "$1" | cut -d' ' -f1)
  printf '%s\t%s\t%s\t%s\t%s\n' \
    "$(date -u +%FT%TZ)" "$h" "$(stat -c%s "$1")" "$2" "${3:-}" >> "$MANIFEST"
}

free_gb(){ df -BG "$1" 2>/dev/null | awk 'NR==2{gsub("G","",$4);print $4}'; }

need_tool(){ command -v "$1" >/dev/null 2>&1; }

# =============================================================================
stage_preflight(){
  c_hdr "PREFLIGHT"
  info "host      : $(hostname)"
  info "cores     : $(nproc)   (using $THREADS; dumps use $DUMP_THREADS)"
  info "memory    : $(free -g | awk '/^Mem:/{print $2}') GB"
  info "project   : $PROJECT"
  info "scratch   : $SCRATCH"
  mkdir -p "$PROJECT" "$SCRATCH"
  local p s; p=$(free_gb "$PROJECT"); s=$(free_gb "$SCRATCH")
  info "free (project): ${p} GB"
  info "free (scratch): ${s} GB"
  echo
  info "Estimated requirements:"
  info "  genomes                 ~12 GB"
  info "  PRJNA342864 .sra+fastq  ~260 GB"
  info "  SRP310563   .sra+fastq  ~340 GB"
  info "  TOTAL                   ~600 GB"
  [[ "${s:-0}" -lt 600 ]] && warn "Under 600 GB free. Use --skip sra-kidney, or --stage by stage."
  echo
  if [[ "$PROJECT" == *" "* ]]; then
    err "Project path contains a SPACE: $PROJECT"
    err "BLAST parses -in/-query as space-separated file LISTS and will break."
    err "Set ACOMYS_PROJECT to a path without spaces."
    return 1
  fi
  ok "no spaces in project path"
  for t in curl wget tar; do
    need_tool "$t" && ok "$t" || warn "$t missing"
  done
  return 0
}

# =============================================================================
stage_setup(){
  c_hdr "SETUP - conda environments"
  if ! need_tool conda; then
    info "installing miniforge to \$HOME/miniforge3 ..."
    wget -q -O /tmp/miniforge.sh \
      "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh" \
      || { err "miniforge download failed"; return 1; }
    bash /tmp/miniforge.sh -b -p "$HOME/miniforge3" >/dev/null || return 1
    # shellcheck disable=SC1091
    source "$HOME/miniforge3/etc/profile.d/conda.sh"
    "$HOME/miniforge3/bin/conda" init bash >/dev/null
    ok "miniforge installed - open a new shell afterwards"
  else
    # shellcheck disable=SC1091
    source "$(conda info --base)/etc/profile.d/conda.sh"
    ok "conda $(conda --version | awk '{print $2}')"
  fi

  if ! conda env list | grep -qE '^acomys-bio\s'; then
    info "creating acomys-bio (blast, sra-tools, hisat2, salmon, hyphy ...)"
    conda create -y -n acomys-bio -c conda-forge -c bioconda \
      ncbi-datasets-cli entrez-direct "sra-tools>=3.0" \
      "blast>=2.15" "diamond>=2.1" "orthofinder>=2.5" \
      "hisat2>=2.2" "stringtie>=2.2" "samtools>=1.18" fastp gffread \
      "salmon>=1.10" fastqc multiqc "hyphy>=2.5" mafft prank \
      >/dev/null 2>&1 && ok "acomys-bio created" || { err "env create failed"; return 1; }
  else
    ok "acomys-bio exists"
  fi

  if ! conda env list | grep -qE '^acomys\s'; then
    info "creating acomys (python analysis stack)"
    conda create -y -n acomys -c conda-forge \
      python=3.11 "numpy<2" pandas scipy scikit-learn statsmodels \
      matplotlib-base seaborn openpyxl requests tqdm biopython pysam \
      >/dev/null 2>&1 && ok "acomys created" || { err "env create failed"; return 1; }
    "$(conda info --base)/envs/acomys/bin/pip" install -q shap dendropy 2>/dev/null
  else
    ok "acomys exists"
  fi

  # Expose bio binaries WITHOUT prepending to PATH - acomys-bio ships its own
  # python and would shadow the analysis interpreter.
  local base abin bbin n=0
  base=$(conda info --base); abin="$base/envs/acomys/bin"; bbin="$base/envs/acomys-bio/bin"
  for t in datasets blastn blastp makeblastdb blastdbcmd diamond orthofinder \
           hyphy mafft prank salmon fastqc multiqc prefetch fasterq-dump \
           vdb-config hisat2 hisat2-build stringtie samtools fastp gffread \
           esearch efetch elink; do
    [[ -x "$bbin/$t" ]] && { ln -sf "$bbin/$t" "$abin/$t"; n=$((n+1)); }
  done
  ok "linked $n tools into the acomys env"

  # SRA toolkit needs a writable cache; point it at scratch.
  mkdir -p "$SCRATCH/ncbi"
  "$bbin/vdb-config" --set "/repository/user/main/public/root=$SCRATCH/ncbi" 2>/dev/null \
    && ok "SRA cache -> $SCRATCH/ncbi" || warn "vdb-config not set; run 'vdb-config -i' if prefetch misbehaves"
  return 0
}

activate(){
  # shellcheck disable=SC1091
  source "$(conda info --base 2>/dev/null || echo "$HOME/miniforge3")/etc/profile.d/conda.sh" 2>/dev/null
  conda activate acomys 2>/dev/null || true
  export PATH="$(conda info --base)/envs/acomys/bin:$PATH"
  export NCBI_API_KEY EMAIL
}

# =============================================================================
stage_genomes(){
  c_hdr "GENOMES - 12 assemblies"
  activate
  need_tool datasets || { err "'datasets' not found - run --stage setup"; return 1; }
  mkdir -p "$REF/genomes"; cd "$REF/genomes" || return 1

  while IFS=$'\t' read -r acc name clade role; do
    [[ -z "$acc" ]] && continue
    if [[ -s "$acc.zip" ]]; then ok "$name  cached ($(du -h "$acc.zip"|cut -f1))"; continue; fi
    info "$name  [$clade] $acc"
    # Unannotated assemblies have no gff3/cds/protein; datasets errors on the
    # whole request rather than skipping, so degrade to genome-only.
    if datasets download genome accession "$acc" \
         --include genome,gff3,rna,cds,protein --filename "$acc.zip" >/dev/null 2>&1; then
      ok "  full package"
    elif datasets download genome accession "$acc" \
         --include genome --filename "$acc.zip" >/dev/null 2>&1; then
      warn "  genome-only (no annotation available)"
    else
      err "  FAILED $acc"; continue
    fi
    record "$acc.zip" "NCBI Datasets $acc" "$name|$clade|$role"
  done <<< "$GENOMES"

  info "unpacking FASTA (needed for BLAST and alignment)"
  for z in *.zip; do
    local a="${z%.zip}"
    [[ -d "$a" ]] && continue
    unzip -qo "$z" -d "$a" 2>/dev/null || python3 -c "
import zipfile,sys; zipfile.ZipFile('$z').extractall('$a')" 2>/dev/null
  done
  ok "genomes: $(ls -1 ./*.zip 2>/dev/null | wc -l)/12   total $(du -sh . | cut -f1)"
  return 0
}

# =============================================================================
stage_geo(){
  c_hdr "GEO - GSE168876 processed tables"
  mkdir -p "$RAW/geo"; cd "$RAW/geo" || return 1
  local base="https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE168876&format=file&file="
  for f in "${GEO_FILES[@]}"; do
    [[ -s "$f" ]] && { ok "$f cached"; continue; }
    local enc="${f//_/%5F}"; enc="${enc//.xlsx/%2Exlsx}"
    curl -fsSL "$base$enc" -o "$f" && { ok "$f"; record "$f" "GEO GSE168876" "processed DE (mouse-ref mapped)"; } \
      || err "failed $f"
  done
  # Two platforms -> one matrix each; there is no combined file.
  for m in "${GEO_MATRICES[@]}"; do
    [[ -s "$m" ]] && { ok "$m cached"; continue; }
    curl -fsSL "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE168nnn/GSE168876/matrix/$m" -o "$m" \
      && { ok "$m"; record "$m" "GEO GSE168876" "series metadata"; } || warn "failed $m"
  done
  return 0
}

# =============================================================================
dump_run(){   # dump_run <SRR> <outdir>
  local run="$1" out="$2"
  mkdir -p "$out" "$SCRATCH/sra_tmp"
  if compgen -G "$out/${run}*.fastq.gz" >/dev/null; then ok "$run FASTQ cached"; return 0; fi

  info "$run: prefetch"
  prefetch "$run" --max-size 200G -O "$SCRATCH/ncbi_dl" >/dev/null 2>&1 \
    || { err "$run prefetch failed"; return 1; }
  local sra; sra=$(find "$SCRATCH/ncbi_dl" -name "$run.sra" | head -1)
  [[ -z "$sra" ]] && { err "$run .sra not found after prefetch"; return 1; }

  info "$run: fasterq-dump ($DUMP_THREADS threads)"
  # -t on fast local disk. Without it fasterq-dump uses the CWD and can die
  # with 'Cannot allocate memory' on network or slow filesystems.
  fasterq-dump "$sra" --split-files --threads "$DUMP_THREADS" --mem 2000MB \
      -t "$SCRATCH/sra_tmp" -O "$out" >/dev/null 2>&1 \
    || { err "$run fasterq-dump failed"; return 1; }

  info "$run: compressing"
  for fq in "$out/$run"*.fastq; do
    [[ -f "$fq" ]] || continue
    if need_tool pigz; then pigz -p "$THREADS" "$fq"; else gzip "$fq"; fi
    record "${fq}.gz" "SRA $run" "raw reads"
  done
  (( KEEP_SRA )) || rm -f "$sra"
  ok "$run done  ($(du -sh "$out" | cut -f1) in $out)"
  return 0
}

stage_sra_utr(){
  c_hdr "SRA - PRJNA342864 (UTR evidence, 2 pooled runs)"
  activate
  info "NOTE: all 15 organs were POOLED before sequencing. This carries NO"
  info "tissue-level expression - it is used ONLY for transcript structure,"
  info "to derive Acomys UTR boundaries without inheriting mouse-projected"
  info "annotation. It is why the spleen arm of this project was dropped."
  for r in "${UTR_RUNS[@]}"; do dump_run "$r" "$SCRATCH/sra/utr"; done
  return 0
}

# =============================================================================
stage_sra_kidney(){
  c_hdr "SRA - SRP310563 (kidney UUO, 18 samples)"
  activate
  info "Needed for WP5: the published tables mapped Acomys reads to the MOUSE"
  info "transcriptome, which under-recovers divergent immune genes. These raw"
  info "reads allow species-native re-quantification."

  local list="$RAW/sra/${KIDNEY_PROJECT}_runs.txt"
  mkdir -p "$RAW/sra"
  if [[ ! -s "$list" ]]; then
    info "resolving run accessions for $KIDNEY_PROJECT ..."
    esearch -db sra -query "$KIDNEY_PROJECT" 2>/dev/null \
      | efetch -format runinfo 2>/dev/null \
      | awk -F, 'NR>1 && $1 ~ /^[SED]RR/ {print $1}' | sort -u > "$list"
  fi
  if [[ ! -s "$list" ]]; then
    err "could not resolve run accessions."
    err "Fetch them manually from:"
    err "  https://www.ncbi.nlm.nih.gov/Traces/study/?acc=$KIDNEY_PROJECT"
    err "then write one SRR per line to: $list"
    return 1
  fi
  local n; n=$(wc -l < "$list")
  ok "$n runs resolved (expect 18)"
  [[ "$n" -ne 18 ]] && warn "expected 18 runs - check $list before trusting it"

  while read -r r; do
    [[ -z "$r" ]] && continue
    dump_run "$r" "$SCRATCH/sra/kidney"
  done < "$list"
  return 0
}

# =============================================================================
summary(){
  c_hdr "SUMMARY"
  info "project : $PROJECT"
  info "scratch : $SCRATCH"
  echo
  [[ -d "$REF/genomes"     ]] && info "genomes      : $(ls -1 "$REF/genomes"/*.zip 2>/dev/null|wc -l)/12  $(du -sh "$REF/genomes" 2>/dev/null|cut -f1)"
  [[ -d "$RAW/geo"         ]] && info "geo tables   : $(ls -1 "$RAW/geo" 2>/dev/null|wc -l) files"
  [[ -d "$SCRATCH/sra/utr" ]] && info "utr reads    : $(ls -1 "$SCRATCH/sra/utr"/*.fastq.gz 2>/dev/null|wc -l) files  $(du -sh "$SCRATCH/sra/utr" 2>/dev/null|cut -f1)"
  [[ -d "$SCRATCH/sra/kidney" ]] && info "kidney reads : $(ls -1 "$SCRATCH/sra/kidney"/*.fastq.gz 2>/dev/null|wc -l) files  $(du -sh "$SCRATCH/sra/kidney" 2>/dev/null|cut -f1)"
  echo
  [[ -f "$MANIFEST" ]] && info "manifest     : $(wc -l < "$MANIFEST") entries -> $MANIFEST"
  info "free now     : project $(free_gb "$PROJECT") GB | scratch $(free_gb "$SCRATCH") GB"
  echo
  info "Add to ~/.bashrc:"
  info "  export ACOMYS_PROJECT=\"$PROJECT\""
  info "  export ACOMYS_SCRATCH=\"$SCRATCH\""
}

# =============================================================================
main(){
  local stages=() skip=()
  while [[ $# -gt 0 ]]; do
    case "$1" in
      --all)       stages=(preflight setup genomes geo sra-utr sra-kidney); shift ;;
      --stage)     stages+=("$2"); shift 2 ;;
      --skip)      skip+=("$2"); shift 2 ;;
      --keep-sra)  KEEP_SRA=1; shift ;;
      --project)   PROJECT="$2"; DATA="$PROJECT/data"; RAW="$DATA/raw"; REF="$DATA/reference"; LOGS="$PROJECT/logs"; MANIFEST="$DATA/fetch_manifest.tsv"; shift 2 ;;
      --scratch)   SCRATCH="$2"; shift 2 ;;
      -h|--help)   sed -n '2,32p' "$0"; exit 0 ;;
      *)           err "unknown option: $1"; exit 1 ;;
    esac
  done
  [[ ${#stages[@]} -eq 0 ]] && { sed -n '2,32p' "$0"; exit 0; }

  mkdir -p "$DATA" "$RAW" "$REF" "$LOGS" "$SCRATCH"
  [[ -f "$MANIFEST" ]] || printf 'utc\tsha256\tbytes\tsource\tnote\n' > "$MANIFEST"
  local logf="$LOGS/fetch_$(date +%Y%m%d_%H%M%S).log"
  info "logging to $logf"

  for s in "${stages[@]}"; do
    [[ " ${skip[*]:-} " == *" $s "* ]] && { warn "skipping $s"; continue; }
    case "$s" in
      preflight)   stage_preflight   ;;
      setup)       stage_setup       ;;
      genomes)     stage_genomes     ;;
      geo)         stage_geo         ;;
      sra-utr)     stage_sra_utr     ;;
      sra-kidney)  stage_sra_kidney  ;;
      *) err "unknown stage: $s" ;;
    esac
    local rc=$?
    (( rc != 0 )) && warn "stage '$s' returned $rc - continuing"
  done 2>&1 | tee -a "$logf"

  summary 2>&1 | tee -a "$logf"
}

main "$@"
