<#
=============================================================================
 ACOMYS PROJECT - NATIVE WINDOWS SETUP
=============================================================================

 Installs every external tool the pipeline needs, natively on Windows.
 NO conda, NO WSL, NO Linux VM - conda is what was crashing the server.

 Run in PowerShell (normal user is fine; admin only helps for Python):

     powershell -ExecutionPolicy Bypass -File setup_windows.ps1

 Options:
     -ToolsDir  C:\Tools      where tools are installed (default C:\Tools)
     -Scratch   D:\scratch    fast local scratch (default C:\acomys_scratch)
     -SkipHeavy               skip SRA Toolkit + kallisto (WP5 only)

 WHAT IT INSTALLS
   BLAST+          WP2 orthology                        official NCBI win64
   ncbi-datasets   WP1 genome download                  official NCBI win64
   SRA Toolkit     WP5 raw reads                        official NCBI win64
   kallisto        WP5 quantification (replaces Salmon) official win release
   PAML codeml     WP4 selection (replaces HyPhy)       official win32
   MAFFT           WP4 alignment                        official windows build

 Python packages are installed with pip into whatever Python is on PATH.
 Install Python 3.11+ from python.org first if you have none.

 NOT INSTALLED, and not needed: HISAT2, StringTie, samtools, Salmon, HyPhy.
 UTR boundaries now come from the curated Acomys dimidiatus annotation
 instead of read alignment - see code/platform_compat.py.
=============================================================================
#>

param(
    [string]$ToolsDir = "C:\Tools",
    [string]$Scratch  = "C:\acomys_scratch",
    [switch]$SkipHeavy
)

$ErrorActionPreference = "Continue"
$ProgressPreference    = "SilentlyContinue"   # much faster downloads

function Hdr($m){ Write-Host "`n=== $m ===" -ForegroundColor Cyan }
function Ok($m) { Write-Host "  [ok] $m"    -ForegroundColor Green }
function Warn($m){Write-Host "  [warn] $m"  -ForegroundColor Yellow }
function Err($m){ Write-Host "  [ERROR] $m" -ForegroundColor Red }
function Info($m){Write-Host "  $m" }

function Get-File($Url, $Dest) {
    if (Test-Path $Dest) { Ok "cached $(Split-Path $Dest -Leaf)"; return $true }
    Info "downloading $(Split-Path $Dest -Leaf)"
    try {
        Invoke-WebRequest -Uri $Url -OutFile $Dest -UseBasicParsing -TimeoutSec 3600
        return $true
    } catch { Err "download failed: $($_.Exception.Message)"; return $false }
}

function Add-ToUserPath($Dir) {
    if (-not (Test-Path $Dir)) { return }
    $cur = [Environment]::GetEnvironmentVariable("Path", "User")
    if ($cur -notlike "*$Dir*") {
        [Environment]::SetEnvironmentVariable("Path", "$cur;$Dir", "User")
        Ok "added to user PATH: $Dir"
    } else { Ok "already on PATH: $Dir" }
    $env:Path = "$env:Path;$Dir"
}

# =============================================================================
Hdr "PREFLIGHT"
Info "tools   : $ToolsDir"
Info "scratch : $Scratch"
New-Item -ItemType Directory -Force -Path $ToolsDir, $Scratch, "$ToolsDir\_dl" | Out-Null

$drive = (Get-Item $Scratch).PSDrive
$freeGB = [math]::Round($drive.Free / 1GB, 0)
Info "free on $($drive.Name): : $freeGB GB"
if ($freeGB -lt 100) { Warn "under 100 GB free - genomes alone need ~30 GB unpacked" }

$py = Get-Command python -ErrorAction SilentlyContinue
if (-not $py) { $py = Get-Command python3 -ErrorAction SilentlyContinue }
if ($py) {
    $ver = & $py.Source --version 2>&1
    Ok "python: $ver  ($($py.Source))"
} else {
    Err "Python not found. Install 3.11+ from https://www.python.org/downloads/"
    Err "Tick 'Add python.exe to PATH' in the installer, then re-run this script."
    exit 1
}

# =============================================================================
Hdr "BLAST+  (WP2 orthology - the critical dependency)"
if (Get-Command blastn -ErrorAction SilentlyContinue) {
    Ok "BLAST+ already installed"
} else {
    $blastUrl = "https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/"
    try {
        $page = Invoke-WebRequest -Uri $blastUrl -UseBasicParsing -TimeoutSec 120
        $file = ($page.Links | Where-Object { $_.href -match 'win64\.tar\.gz$' } |
                 Select-Object -First 1).href
    } catch { $file = $null }
    if ($file) {
        $dl = "$ToolsDir\_dl\$file"
        if (Get-File "$blastUrl$file" $dl) {
            Info "extracting BLAST+"
            tar -xzf $dl -C $ToolsDir 2>$null
            $bin = Get-ChildItem $ToolsDir -Directory -Filter "ncbi-blast-*" |
                   Select-Object -First 1
            if ($bin) { Add-ToUserPath "$($bin.FullName)\bin"; Ok "BLAST+ installed" }
            else { Err "extraction produced no ncbi-blast-* folder" }
        }
    } else {
        Err "could not resolve the BLAST+ win64 archive."
        Err "Download manually: https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/"
    }
}

# =============================================================================
Hdr "NCBI datasets  (WP1 genome download)"
if (Get-Command datasets -ErrorAction SilentlyContinue) {
    Ok "datasets already installed"
} else {
    $d = "$ToolsDir\ncbi-datasets"
    New-Item -ItemType Directory -Force -Path $d | Out-Null
    $u = "https://ftp.ncbi.nlm.nih.gov/pub/datasets/command-line/v2/windows-amd64/datasets.exe"
    if (Get-File $u "$d\datasets.exe") { Add-ToUserPath $d; Ok "datasets installed" }
}

# =============================================================================
Hdr "PAML  (WP4 selection tests - replaces HyPhy)"
if (Get-Command codeml -ErrorAction SilentlyContinue) {
    Ok "codeml already installed"
} else {
    $u = "https://github.com/abacus-gene/paml/archive/refs/heads/master.zip"
    $dl = "$ToolsDir\_dl\paml.zip"
    if (Get-File $u $dl) {
        Expand-Archive -Path $dl -DestinationPath $ToolsDir -Force
        $pb = Get-ChildItem $ToolsDir -Directory -Filter "paml*" | Select-Object -First 1
        if ($pb -and (Test-Path "$($pb.FullName)\bin")) {
            Add-ToUserPath "$($pb.FullName)\bin"; Ok "PAML installed"
        } else { Warn "PAML extracted but bin\ not found - check $ToolsDir" }
    }
}

# =============================================================================
Hdr "MAFFT  (WP4 alignment)"
if (Get-Command mafft -ErrorAction SilentlyContinue) {
    Ok "mafft already installed"
} else {
    Info "MAFFT ships as an installer; fetch from:"
    Info "  https://mafft.cbrc.jp/alignment/software/windows.html"
    Info "  (all-in-one package, then add its folder to PATH)"
    Warn "skipped - only needed for WP4"
}

# =============================================================================
if (-not $SkipHeavy) {
    Hdr "SRA Toolkit  (WP5 raw reads)"
    if (Get-Command fasterq-dump -ErrorAction SilentlyContinue) {
        Ok "SRA Toolkit already installed"
    } else {
        $u = "https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/current/sratoolkit.current-win64.zip"
        $dl = "$ToolsDir\_dl\sratoolkit.zip"
        if (Get-File $u $dl) {
            Expand-Archive -Path $dl -DestinationPath $ToolsDir -Force
            $sb = Get-ChildItem $ToolsDir -Directory -Filter "sratoolkit*" |
                  Select-Object -First 1
            if ($sb) { Add-ToUserPath "$($sb.FullName)\bin"; Ok "SRA Toolkit installed" }
        }
    }

    Hdr "kallisto  (WP5 quantification - replaces Salmon)"
    if (Get-Command kallisto -ErrorAction SilentlyContinue) {
        Ok "kallisto already installed"
    } else {
        Info "kallisto Windows builds are on the releases page:"
        Info "  https://github.com/pachterlab/kallisto/releases"
        Info "  download the windows zip, extract to $ToolsDir\kallisto, add to PATH"
        Warn "skipped - only needed for WP5"
    }
}

# =============================================================================
Hdr "PYTHON PACKAGES"
$pkgs = @("numpy<2","pandas","scipy","scikit-learn","statsmodels","matplotlib",
          "seaborn","biopython","openpyxl","requests","tqdm","shap","dendropy")
Info "installing: $($pkgs -join ', ')"
& $py.Source -m pip install --quiet --upgrade pip 2>&1 | Out-Null
& $py.Source -m pip install --quiet $pkgs 2>&1 | Out-Null
if ($LASTEXITCODE -eq 0) { Ok "python packages installed" }
else { Warn "some packages may have failed - check with: python -m pip list" }

# pydeseq2 replaces R limma/edgeR for WP5; optional, can be slow to build
& $py.Source -m pip install --quiet pydeseq2 2>&1 | Out-Null

# =============================================================================
Hdr "ENVIRONMENT"
[Environment]::SetEnvironmentVariable("ACOMYS_SCRATCH", $Scratch, "User")
$env:ACOMYS_SCRATCH = $Scratch
Ok "ACOMYS_SCRATCH = $Scratch"

Hdr "VERIFY"
Info "Open a NEW PowerShell window (so PATH refreshes), then:"
Info ""
Info "  cd '$PSScriptRoot'"
Info "  python code\platform_compat.py"
Info "  python run_all.py --check"
Info ""
Info "Then fetch data and run the primary analysis:"
Info "  python code\01_fetch_data.py --step genomes"
Info "  python code\02_build_ortholog_scaffold.py --all"
Info "  python code\03_kmer_phylo_controlled.py --arm residual --region promoter"
Info ""
Write-Host "=== setup complete ===" -ForegroundColor Cyan
