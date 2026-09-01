#!/usr/bin/env python3
"""
================================================================================
SHARED CONFIGURATION - single source of truth
Acomys cahirinus regeneration project
================================================================================

Every other script imports from here. Nothing is duplicated. If an accession,
gene set, seed or threshold needs to change, it changes once, in this file.

Import from any pipeline script with a plain:

    import config as cfg

(Deliberately NOT named 00_config.py -- a leading digit makes normal import
illegal and forces an importlib workaround that breaks dataclasses. Numbered
filenames are for pipeline STEPS; config is not a step.)

VERIFICATION STATUS: all assembly accessions below were confirmed against NCBI
on 13 August 2026 by the WP1 species annotation audit. See
WP1_species_annotation_audit.md for the record.
================================================================================
"""

from __future__ import annotations

import random
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

# ==============================================================================
# REPRODUCIBILITY
# ==============================================================================

RANDOM_SEED = 42


def set_all_seeds(seed: int = RANDOM_SEED):
    """Call at the top of every script. Do not set seeds anywhere else."""
    import os
    random.seed(seed)
    np.random.seed(seed)
    os.environ["PYTHONHASHSEED"] = str(seed)


# ==============================================================================
# PATHS
# ==============================================================================

PROJECT_ROOT = Path(__file__).resolve().parent.parent
CODE = PROJECT_ROOT / "code"
DATA = PROJECT_ROOT / "data"
RAW = DATA / "raw"
REF = DATA / "reference"
ORTHO = REF / "orthologs"
RESULTS = PROJECT_ROOT / "results"
LOGS = PROJECT_ROOT / "logs"
MANIFEST = DATA / "provenance_manifest.json"

ALL_DIRS = [RAW / "geo", RAW / "sra", REF / "genomes", REF / "transcriptome",
            ORTHO, RESULTS, LOGS]


def ensure_dirs():
    for d in ALL_DIRS:
        d.mkdir(parents=True, exist_ok=True)


# ==============================================================================
# NCBI
# ==============================================================================

ENTREZ_EMAIL = "ahmad.hassanat@gmail.com"


def load_dotenv(path: Path | None = None) -> dict:
    """Read KEY=value pairs from a gitignored .env into os.environ.

    Credentials live in .env, never in this file. config.py is committed and
    the project is destined for a public repository; a key pasted in here goes
    public with it and cannot be un-published. .env is covered by .gitignore.

    Existing environment variables win, so an exported value overrides .env.
    """
    import os
    path = path or (PROJECT_ROOT / ".env")
    loaded = {}
    if not path.exists():
        return loaded
    for line in path.read_text().splitlines():
        line = line.strip()
        if not line or line.startswith("#") or "=" not in line:
            continue
        k, v = line.split("=", 1)
        k, v = k.strip(), v.strip().strip('"').strip("'")
        loaded[k] = v
        os.environ.setdefault(k, v)
    return loaded


load_dotenv()


def ncbi_api_key() -> str | None:
    """The NCBI key, or None. Validated - a malformed key makes EVERY Entrez
    request return HTTP 400, which reads as missing data rather than a bad
    credential."""
    import os
    key = os.environ.get("NCBI_API_KEY")
    if not key:
        return None
    key = key.strip()
    if len(key) < 30 or not all(c in "0123456789abcdefABCDEF" for c in key):
        log(f"NCBI_API_KEY looks malformed (len {len(key)}, expected ~36 hex "
            f"chars). Ignoring it - an invalid key makes every Entrez call "
            f"fail with HTTP 400.", "WARN")
        os.environ.pop("NCBI_API_KEY", None)
        return None
    return key


# ==============================================================================
# SPECIES PANEL - VERIFIED 13 Aug 2026
# ==============================================================================
# Phylogeny (Steppan et al. 2004; Alhajeri et al. 2015; BMC Biology 2024):
#
#        +-- Deomyinae     Acomys                      REGENERATIVE
#   +----+
#   |    +-- Gerbillinae   Meriones                    SISTER CLADE
# --+
#   +-- Murinae            Mus, Rattus, Apodemus, ...  non-regenerative
#   |
#   +-- Cricetidae         Peromyscus, Mesocricetus    outgroup / root
#
# Deomyinae + Gerbillinae split from Murinae ~18.3 Ma.

@dataclass
class Species:
    name: str
    taxid: int
    clade: str
    regenerative: str            # "yes" | "no" | "unknown"
    assembly: str | None
    assembly_name: str | None
    annotation: str              # "refseq" | "projected" | "none"
    divergence_ma: float         # from Acomys cahirinus
    usable: bool
    note: str = ""


SPECIES_PANEL: list[Species] = [
    # ---- Deomyinae (focal clade) ----
    Species("Acomys cahirinus", 10068, "Deomyinae", "yes",
            "GCA_029890205.1", "ASM2989020v1", "projected", 0.0, True,
            "FOCAL. Nanopore+Hi-C, N50 ~127 Mb, ~98.5% complete. "
            "Annotation is CAT-projected from GRCm39 - see circularity caveat."),
    Species("Acomys russatus", 60746, "Deomyinae", "yes",
            "GCF_903995435.1", "mAcoRus1.1", "refseq", 5.0, True,
            "RefSeq-annotated. Use as INDEPENDENT check on A. cahirinus "
            "UTR boundaries - it does not inherit mouse-projected structure."),
    Species("Acomys dimidiatus", 269703, "Deomyinae", "yes",
            "GCA_907164435.1", "mAcoDim1_REL_1905", "ensembl", 2.0, True,
            "CORRECTED 13 Aug 2026: annotation DOES exist, via Ensembl "
            "(manually curated), not NCBI RefSeq. `datasets` will not supply "
            "a GFF - WP2 must fetch gene models from Ensembl or project them "
            "from A. russatus."),

    # ---- Gerbillinae (sister clade) ----
    Species("Meriones unguiculatus", 10047, "Gerbillinae", "no",
            "GCF_030254825.1", "Bangor_MerUng_6.1", "refseq", 15.0, True,
            "Nearest annotated non-regenerator."),
    Species("Psammomys obesus", 48139, "Gerbillinae", "no",
            "GCF_907164565.1", "mPsaObe1.curated_primary_1811", "refseq",
            15.0, True,
            "CORRECTED 13 Aug 2026: RefSeq-annotated chromosome-level Sanger "
            "assembly DOES exist. The original audit searched the literature "
            "and found only a mitochondrial genome. This resolves the former "
            "'Gerbillinae n=1' design gap."),

    # ---- Murinae ----
    Species("Mus musculus", 10090, "Murinae", "no",
            "GCF_000001635.27", "GRCm39", "refseq", 18.3, True,
            "Standard comparator; supplies the injury data (GSE168876)."),
    Species("Rattus norvegicus", 10116, "Murinae", "no",
            "GCF_015227675.2", "mRatBN7.2", "refseq", 18.3, True, ""),
    Species("Apodemus sylvaticus", 10129, "Murinae", "no",
            "GCF_947179515.1", "mApoSyl1.1", "refseq", 18.3, True,
            "Darwin Tree of Life, chromosome-level."),
    Species("Grammomys surdaster", 491861, "Murinae", "no",
            "GCF_004785775.1", "NIH_TR_1.0", "refseq", 18.3, True,
            "RefSeq AR100, 21,298 coding genes. Status resolved from "
            "'unknown' to 'no': Riddell et al. 2025 PNAS found ALL "
            "non-deomyine rodents heal by fibrotic repair."),
    Species("Mastomys coucha", 35658, "Murinae", "no",
            "GCF_008632895.1", "UCSF_Mcou_1", "refseq", 18.3, True,
            "RefSeq AR100, 20,153 coding genes. Status resolved as above. "
            "NOTE: `datasets` served this one genome-only - no annotation in "
            "the package despite the RefSeq listing. WP2 must handle it."),

    # ---- Cricetidae (outgroup) ----
    Species("Peromyscus maniculatus", 10042, "Cricetidae", "no",
            "GCF_049852395.1", "HU_Pman_BW_mat_3.1", "refseq", 25.0, True,
            "Roots the tree."),
    Species("Mesocricetus auratus", 10036, "Cricetidae", "no",
            "GCF_017639785.1", "BCM_Maur_2.0", "refseq", 25.0, True,
            "Roots the tree."),
]

USABLE_SPECIES = [s for s in SPECIES_PANEL if s.usable]
FOCAL_SPECIES = "Acomys cahirinus"
MIN_SPECIES_REQUIRED = 8

DIVERGENCE_MA = {s.name: s.divergence_ma for s in SPECIES_PANEL}
CLADE = {s.name: s.clade for s in SPECIES_PANEL}
REGENERATIVE = {s.name: (s.regenerative == "yes") for s in SPECIES_PANEL}
ASSEMBLY = {s.name: s.assembly for s in SPECIES_PANEL if s.assembly}

# ------------------------------------------------------------------------------
# KNOWN DESIGN GAPS - carry these into the manuscript limitations, do not
# quietly forget them.
# ------------------------------------------------------------------------------
DESIGN_GAPS = {
    "regeneration_is_a_clade_trait": (
        "Riddell et al. 2025 (PNAS 122:e2420726122) tested ear-pinna "
        "regeneration across 11 rodents and found that ONLY the deomyines - "
        "Acomys spp. and Lophuromys zena - rebuild complex tissue; every "
        "non-deomyine heals by fibrotic repair. Regeneration is therefore a "
        "CLADE trait of Deomyinae, not an Acomys peculiarity. Two "
        "consequences: (i) the regenerative/non-regenerative labels in this "
        "panel now rest on published phenotyping rather than assumption; "
        "(ii) all three regenerative species here are Acomys, so "
        "'regenerative' remains confounded with 'genus Acomys'. Lophuromys "
        "would break that confound but has no nuclear genome assembly. Flag "
        "this as the single most valuable species to add if one appears."
    ),
    "acomys_annotation_circularity": (
        "A. cahirinus annotation is CAT-projected FROM GRCm39. Using "
        "mouse-projected gene models to measure divergence FROM mouse is "
        "circular, and would bias UTR boundaries toward mouse-like structure "
        "- precisely the measurement WP3 makes. MITIGATION: (i) validate "
        "A. cahirinus UTR boundaries against A. russatus (GCF_903995435.1), "
        "which is RefSeq-annotated independently; (ii) define UTRs de novo "
        "from PRJNA342864 transcript evidence; (iii) report WP3 results "
        "under all three annotation sources."
    ),
    "clade_imbalance": (
        "Usable species per clade after the 13 Aug 2026 corrections: "
        "Deomyinae 3, Gerbillinae 2, Murinae 5, Cricetidae 2. Murinae is "
        "still over-represented. Weight the PGLS accordingly and report "
        "leave-one-species-out sensitivity for every species."
    ),
    "no_spleen_data_at_all": (
        "SCOPE CHANGE 13 Aug 2026. There is NO Acomys spleen expression data "
        "in any public dataset. PRJNA342864 pooled all 15 organs before "
        "sequencing (2 runs), so it has no tissue resolution; GSE168876 is "
        "kidney only. The spleen arm was dropped and the project is now "
        "kidney-only. Myeloid involvement is addressed by deconvolving the "
        "kidney RNA-seq, which infers the COMPOSITION of the infiltrate but "
        "says nothing about where those cells came from. Do not write "
        "'spleen', 'reservoir' or 'trafficking' anywhere in the manuscript."
    ),
}


# ==============================================================================
# GENE SETS
# ==============================================================================

FOCAL_GENES = {
    "CX3CL1": "Fractalkine. Membrane-tethered; recruits/retains patrolling monocytes.",
    "CX3CR1": "Fractalkine receptor. Marks resident/patrolling, non-fibrogenic macrophages.",
    "CCR2":   "CCL2 receptor. Ly6C-hi monocyte recruitment; pro-fibrotic axis.",
    "TGFB1":  "Master pro-fibrotic cytokine. Drives myofibroblast transition.",
    "TGFB3":  "Anti-scarring isoform. Elevated in scarless fetal wound healing.",
}

# ------------------------------------------------------------------------------
# EXCLUDED FROM THE FOCAL MODULE - decided 14 Aug 2026, BEFORE any WP3 analysis
# ------------------------------------------------------------------------------
EXCLUDED_GENES = {
    "CCL2": (
        "MCP-1. Dropped after WP2 orthology QC, before any analysis was run. "
        "Evidence: (i) no annotated ortholog in the congener A. russatus, so "
        "Tier 2 had no query; (ii) recovered CDS in A. cahirinus was 112 bp "
        "against ~447 expected; (iii) implausible promoter lengths in four "
        "species (213-6557 bp for a 2000 bp query); (iv) native symbol lookup "
        "failed in most of the panel even with aliases; (v) CCL2 lies in a "
        "rapidly evolving chemokine cluster with lineage-specific "
        "duplications, so a homology hit is as likely to be a paralog as the "
        "ortholog. A wrong ortholog is worse than a missing one because it "
        "looks fine. The CCL2-CCR2 recruitment axis is still represented by "
        "CCR2, which recovered at 1.00x across CDS, 3'UTR and promoter."
    ),
}

CONTROL_SETS = {
    "housekeeping": ["ACTB", "GAPDH", "TUBB5", "RPL13A", "PPIA", "B2M",
                     "SDHA", "TBP", "HPRT1", "YWHAZ"],
    "immune_nonfibrotic": ["CD19", "MS4A1", "CD3E", "IL7R", "TLR4", "TLR9",
                           "NLRP3", "IFNG", "IL2", "CD8A"],
    "fibrosis_effector": ["COL1A1", "COL3A1", "ACTA2", "FN1", "TIMP1",
                          "MMP2", "MMP9", "SERPINE1", "CTGF", "SMAD3"],
    "random_background": [],   # 200 genes drawn by 02_build_ortholog_scaffold.py
}

N_RANDOM_BACKGROUND = 200

# ------------------------------------------------------------------------------
# EXTENDED PANEL - secondary endpoint (E. Hassanat, received 20 Aug 2026 22:50)
# ------------------------------------------------------------------------------
# An a priori selection from canonical renal injury/fibrosis literature
# (exemplar: Lan HY 2011, Int J Biol Sci 7(7):1056-1067), NOT a commercial kit.
# The selector attests she had not seen the RNA-seq data. Logged verbatim in
# PREREGISTRATION.md before any analysis was run against it, and BEFORE any
# discussion of testing design, so that the design discussion could not shape
# panel membership.
#
# DO NOT add, remove, substitute or reorder members of this panel. If a gene
# fails to resolve, that is reported as a result - it is not repaired.
#
# The three blocks are analysed SEPARATELY (E.H. decision, 20 Aug). They form a
# mechanistic gradient - upstream signalling -> injury/repair readout ->
# downstream effectors - so testing 5+5+5 locates the boundary of the effect,
# which a single 15-gene average would blur.
EXTENDED_PANEL = {
    "ECM_fibrosis":          ["COL1A1", "COL3A1", "FN1", "TIMP1", "ACTA2"],
    "inflammatory_axis":     ["CCL2", "CXCL2", "ICAM1", "TNF", "IL6"],
    "tubular_damage_repair": ["HAVCR1", "LCN2", "SPP1", "MYC", "EGF"],
}

# Reported alongside the blocks in supplementary tables (E.H. decision).
EXTENDED_15 = [g for v in EXTENDED_PANEL.values() for g in v]
PANEL_20 = list(FOCAL_GENES) + EXTENDED_15

# Known overlap with sets already under test. Recorded so the manuscript
# methods can state it, and so H1/H2 can exclude the affected control set.
# All five ECM genes are members of the fibrosis_effector CONTROL set, so a
# Fisher-style ECM-vs-fibrosis_effector comparison would pit a set against one
# sharing half its members. E.H. decision (20 Aug): test the ECM block against
# housekeeping and immune_nonfibrotic only, plus the all-gene background, and
# document the overlap.
PANEL_CONTROL_OVERLAP = {
    "ECM_fibrosis": ["fibrosis_effector"],
    "inflammatory_axis": ["macro_proinflam"],   # CXCL2, TNF, IL6
    "tubular_damage_repair": [],                # no overlap
}


def controls_for(block: str) -> list[str]:
    """Control sets usable for a given panel block, overlap removed."""
    bad = set(PANEL_CONTROL_OVERLAP.get(block, []))
    return [k for k, v in CONTROL_SETS.items() if v and k not in bad]

# ------------------------------------------------------------------------------
# SYMBOL ALIASES
# ------------------------------------------------------------------------------
# NCBI Gene symbols are not stable across rodent species. CCL2 in particular
# appears under several historical names, which is why a symbol-only lookup
# failed in 9 of 12 species on the first WP2 run. Aliases are tried in order
# after the primary symbol.
#
# Verify any alias-derived hit by synteny before trusting it: an alias can
# match a paralog in a species where the true ortholog is absent.
GENE_ALIASES = {
    "CCL2":   ["Ccl2", "Mcp1", "MCP-1", "Scya2", "Je", "Sigje", "MCAF"],
    "CCR2":   ["Ccr2", "Cmkbr2", "Ck-r2", "Cd192"],
    "CX3CL1": ["Cx3cl1", "Scyd1", "Fractalkine", "Neurotactin", "Abcd3"],
    "CX3CR1": ["Cx3cr1", "Ccrl1", "V28", "Gpr13", "Cmkbrl1"],
    "TGFB1":  ["Tgfb1", "Tgfb-1", "TGF-beta1"],
    "TGFB3":  ["Tgfb3", "Tgfb-3", "TGF-beta3"],
    "CTGF":   ["Ctgf", "Ccn2", "Fisp12", "Hcs24"],
    "SERPINE1": ["Serpine1", "Pai1", "Planh1", "Pai-1"],
    "TUBB5":  ["Tubb5", "Tubb", "M-beta-5"],
    "ACTA2":  ["Acta2", "Actsa", "alpha-SMA"],
}


def symbol_candidates(gene: str) -> list[str]:
    """Primary symbol first, then documented aliases."""
    out = [gene]
    for a in GENE_ALIASES.get(gene, []):
        if a.upper() != gene.upper():
            out.append(a)
    return out


# Genes needing manual orthology verification before use.
ORTHOLOGY_WATCHLIST = {
    "CCL2": "Sits in a rapidly evolving chemokine cluster with lineage-specific "
            "duplications. Verify synteny manually; exclude if unresolvable.",
    "CX3CL1": "Check for tandem paralogs in non-model species.",
    "TGFB3": "Confirm not confused with TGFB2 by automated assignment.",
}


def all_genes() -> list[str]:
    out = list(FOCAL_GENES)
    for genes in CONTROL_SETS.values():
        out.extend(genes)
    return out


# ==============================================================================
# SEQUENCE REGIONS
# ==============================================================================

REGIONS = ["promoter", "utr5", "cds", "utr3", "full"]
PROMOTER_UPSTREAM_BP = 2000
MIN_REGION_LENGTH = 60          # shorter regions are dropped as stubs
FALLBACK_UTR_LENGTH = 500       # window from CDS boundary where UTR unannotated


# ==============================================================================
# ANALYSIS PARAMETERS
# ==============================================================================

K_SWEEP = [3, 4, 5, 6]
DEFAULT_K = 5
N_PERMUTATIONS = 1000
N_BOOTSTRAP = 2000
N_NULL_SHUFFLES = 200
CV_FOLDS = 5
ALPHA = 0.05

# Pre-specified decision rule for H1. Do not relax this after seeing results.
H1_DECISION_RULE = (
    "Claim divergence ONLY IF the focal module exceeds ALL THREE control sets "
    "at q < 0.05 (BH within arm) AND survives the Arm 2 composition null "
    "(p_empirical < 0.05)."
)


# ==============================================================================
# DATASETS
# ==============================================================================

GSE = "GSE168876"
GSE_BIOPROJECT = "PRJNA714406"
GSE_SRA = "SRP310563"

# Two platforms -> GEO splits the series matrix into per-platform files.
# There is no single GSE168876_series_matrix.txt.gz.
GSE_PLATFORMS = {
    "GPL24247": "Illumina NovaSeq 6000 (Mus musculus)",
    "GPL29848": "Illumina NovaSeq 6000 (Acomys cahirinus)",
}

GSE168876_SAMPLES = [
    ("GSM5171802", "Mus musculus",     "sham", 0), ("GSM5171803", "Mus musculus", "sham", 0),
    ("GSM5171804", "Mus musculus",     "sham", 0),
    ("GSM5171805", "Mus musculus",     "UUO",  2), ("GSM5171806", "Mus musculus", "UUO", 2),
    ("GSM5171807", "Mus musculus",     "UUO",  2),
    ("GSM5171808", "Mus musculus",     "UUO",  5), ("GSM5171809", "Mus musculus", "UUO", 5),
    ("GSM5171810", "Mus musculus",     "UUO",  5),
    ("GSM5171811", "Acomys cahirinus", "sham", 0), ("GSM5171812", "Acomys cahirinus", "sham", 0),
    ("GSM5171813", "Acomys cahirinus", "sham", 0),
    ("GSM5171814", "Acomys cahirinus", "UUO",  2), ("GSM5171815", "Acomys cahirinus", "UUO", 2),
    ("GSM5171816", "Acomys cahirinus", "UUO",  2),
    ("GSM5171817", "Acomys cahirinus", "UUO",  5), ("GSM5171818", "Acomys cahirinus", "UUO", 5),
    ("GSM5171819", "Acomys cahirinus", "UUO",  5),
]

# ------------------------------------------------------------------------------
# PRJNA342864 - Mamrot et al. 2017, Sci Rep. "15-organ transcriptome".
#
# CORRECTED 13 Aug 2026 after reading the methods. The 15 organs were POOLED
# BEFORE SEQUENCING ("were pooled for sequencing... not multiplexed at the time
# of sequencing due to cost"), yielding just TWO SRA runs - male and female
# pools. The 2023 genome paper independently describes it as "bulk RNA-seq from
# multiple pooled organs".
#
# CONSEQUENCE: there is NO per-organ Acomys expression data. No spleen
# expression, baseline or otherwise. The spleen arm was dropped entirely on
# 13 Aug 2026; this project is now kidney-only. See PREREGISTRATION.md sec.9.
#
# STILL USEFUL FOR: transcript STRUCTURE. Pooled tissue is fine for calling UTR
# boundaries from read coverage, which is what WP2 needs to escape the
# mouse-projected annotation circularity. Expression is not required for that.
# ------------------------------------------------------------------------------
ACOMYS_RNASEQ_BIOPROJECT = "PRJNA342864"
ACOMYS_RNASEQ_RUNS = ["SRR4279903", "SRR4279904"]   # male pool, female pool
ACOMYS_RNASEQ_NOTE = (
    "Pooled 15-organ RNA-seq. Use for UTR structure only - it carries NO "
    "tissue-level expression information."
)

# Other Acomys RNA-seq BioProjects noted in the literature. Unverified; listed
# as leads for WP6 cross-tissue replication, not relied upon.
# ==============================================================================
# SINGLE-CELL MYELOID DATA - found 18 Aug 2026 by direct NCBI query
# ==============================================================================
# GSE182141 (Adam, Potter & Seifert; PMID 38228141, 39990382)
#   "A subset of pro-inflammatory spiny mice macrophages promotes a
#    regenerative phenotype in mouse fibroblasts"
#
#   10X single-cell RNA-seq, BOTH species, injured EAR PINNAE, timecourse
#   d0/3/5/10/15. Processed CellRanger matrices (903 MB) - no alignment needed.
#   BioProject PRJNA754776 | SRA SRP332607
#
# WHY IT MATTERS: H4 proposed INFERRING myeloid composition by deconvolving
# n=3 bulk kidney samples. This OBSERVES it at single-cell resolution across
# five timepoints in both species. It is a far stronger test of the same
# hypothesis.
#
# THE CAVEAT: ear pinnae, not kidney. Different organ and injury model. The
# ear punch is the canonical assay defining the deomyine regenerative
# phenotype (Riddell et al. 2025), so it is defensible - but H4 would then be
# answered in ear, and the manuscript must say so plainly rather than eliding
# the difference.
SINGLE_CELL_MYELOID = {
    "GSE182141": {
        "bioproject": "PRJNA754776",
        "sra": "SRP332607",
        "tissue": "ear pinnae (NOT kidney)",
        "species": ["Acomys cahirinus", "Mus musculus"],
        "design": "10X scRNA-seq, d0/3/5/10/15 post-injury, both species",
        "samples": ["GSM5519169 Mus00", "GSM5519170 Mus03", "GSM5519171 Mus05",
                    "GSM5519172 Mus10", "GSM5519173 Mus15",
                    "GSM5519174 Acomys00", "GSM5519175 Acomys03",
                    "GSM5519176 Acomys05", "GSM5519177 Acomys10",
                    "GSM5519178 Acomys15"],
        "processed": "GSE182141_RAW.tar (903.6 MB, MTX/TSV CellRanger)",
        "url": "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE182141&format=file",
    },
    "GSE224879": {"note": "ERK-dependent switch antagonising fibrosis; scRNA-seq n=17"},
    "GSE224433": {"note": "companion to GSE224879; scRNA-seq n=6"},
    "GSE216723": {"note": "injured Mus and Acomys ear d0, 10X"},
}

# Additional Acomys assemblies found 18 Aug 2026. Deomyinae could go 3 -> 6.
ACOMYS_ADDITIONAL_ASSEMBLIES = {
    "GCF_964271855.1": "Acomys minous mAcoMin1.4 - REFSEQ ANNOTATED. A second "
                       "independently annotated congener; could supply the "
                       "full-length 3'UTRs that H1 lacked.",
    "GCA_964186705.1": "Acomys percivali mAcoPer2-scaffold",
    "GCA_964186715.1": "Acomys kempi mAcoKem2_scaffolds",
    "GCA_030555615.1": "Acomys cahirinus ASM3055561v1 - second assembly of the "
                       "focal species; cross-check for the CAT annotation.",
}

ACOMYS_RNASEQ_LEADS = {
    "PRJNA184055": "unverified - speciation / microclimate context",
    "PRJNA292021": "unverified - lifespan / skin regeneration context",
}

# ==============================================================================
# KEY LITERATURE - engage with these in the manuscript
# ==============================================================================

KEY_PAPERS = {
    "Riddell et al. 2025 PNAS 122:e2420726122": (
        "Regeneration across 11 rodents. Only deomyines (Acomys, Lophuromys) "
        "regenerate; all non-deomyines scar. FOUNDATIONAL for this project's "
        "phylogenetic framing and the source of our phenotype labels."
    ),
    "Okamura & Majesky 2021 iScience (GSE168876)": (
        "Source of the kidney injury data. Concluded responses were LARGELY "
        "CONSERVED between species - must be engaged with directly, not "
        "contradicted naively."
    ),
    "Nguyen et al. 2023 G3 jkad177": (
        "Acomys chromosome-scale assembly GCA_029890205.1. Note the CAT "
        "projection from GRCm39 - the circularity problem."
    ),
    "Mamrot et al. 2017 Sci Rep 7:8996": (
        "PRJNA342864. Organs POOLED before sequencing - no tissue resolution."
    ),
    "TO READ - Science 2025, 'Reactivation of mammalian regeneration by "
    "turning on an evolutionarily disabled genetic switch'": (
        "doi 10.1126/science.adp0176. Potentially highly relevant to the "
        "regulatory-divergence hypothesis. Not yet assessed."
    ),
    "TO READ - bioRxiv 2025.02.10.637521": (
        "'Specific cell states underlie complex tissue regeneration in spiny "
        "mice'. Single-cell; may bear on the myeloid deconvolution arm."
    ),
}

# The caveat that drives WP5-C1.
GSE_CAVEAT = (
    "Authors quantified Acomys reads against the MOUSE transcriptome. "
    "Cross-species pseudo-alignment under-recovers divergent genes, which are "
    "disproportionately chemokines and immune genes. Published DE tables are "
    "ORIENTATION ONLY. All divergence claims must rest on species-native "
    "re-quantification."
)


# ==============================================================================
# PLOTTING
# ==============================================================================

PLOT_STYLE = {
    "figure.dpi": 150, "savefig.dpi": 300,
    "font.family": "sans-serif", "font.size": 9,
    "axes.spines.top": False, "axes.spines.right": False,
}
COLOR_FOCAL = "#c0392b"
COLOR_CONTROL = "#95a5a6"
CLADE_COLORS = {"Deomyinae": "#c0392b", "Gerbillinae": "#e67e22",
                "Murinae": "#2980b9", "Cricetidae": "#7f8c8d"}


def log(msg: str, level: str = "INFO"):
    print(f"[{level}] {msg}", flush=True)


def summary():
    """Print the panel status. Run this to re-check the WP2 gate."""
    log("=" * 74)
    log("SPECIES PANEL STATUS")
    log("=" * 74)
    by_clade: dict[str, list[Species]] = {}
    for s in SPECIES_PANEL:
        by_clade.setdefault(s.clade, []).append(s)
    for clade, members in by_clade.items():
        n_ok = sum(1 for m in members if m.usable)
        log(f"{clade:14s} {n_ok}/{len(members)} usable")
        for m in members:
            mark = "OK  " if m.usable else "SKIP"
            log(f"   [{mark}] {m.name:24s} {m.assembly or '-':20s} {m.annotation}")
    n = len(USABLE_SPECIES)
    log("-" * 74)
    log(f"TOTAL USABLE: {n} / {len(SPECIES_PANEL)}   (minimum required: {MIN_SPECIES_REQUIRED})")
    log(f"GATE: {'PASS - PGLS design viable' if n >= MIN_SPECIES_REQUIRED else 'FAIL - invoke fallback'}")
    log("-" * 74)
    log("DESIGN GAPS carried forward:")
    for key, text in DESIGN_GAPS.items():
        log(f"  * {key}")
        for line in _wrap(text, 68):
            log(f"      {line}")
    log("=" * 74)


def _wrap(text: str, width: int) -> list[str]:
    words, lines, cur = text.split(), [], ""
    for w in words:
        if len(cur) + len(w) + 1 > width:
            lines.append(cur)
            cur = w
        else:
            cur = f"{cur} {w}".strip()
    if cur:
        lines.append(cur)
    return lines


if __name__ == "__main__":
    summary()
