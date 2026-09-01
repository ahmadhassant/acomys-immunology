# acomys-immunology — analysis code

Analysis pipeline for a pre-registered, phylogeny-controlled comparative study
of the CX3CL1–CX3CR1–CCR2–TGF-β module across murid rodents, using public data
only. Fully computational; no wet-lab component.

> **Status: code released ahead of the manuscript.**
> This repository currently contains the analysis code and environment
> specifications only. The pre-registration, result tables, figures and the
> manuscript are withheld until the paper is submitted, at which point they
> will be added here and archived with a DOI.

---

## Quick start

Linux or WSL2. bioconda has no Windows builds, so a native Windows conda
environment cannot run the alignment, selection or quantification steps.

```bash
# Two environments, deliberately separate
conda env create -f environment.yml        # acomys      - Python stack
conda env create -f environment-bio.yml    # acomys-bio  - CLI binaries
conda activate acomys

# Expose the binaries WITHOUT exposing acomys-bio's interpreter
for t in datasets blastn blastp makeblastdb diamond orthofinder \
         hyphy mafft prank salmon fastqc multiqc prefetch fasterq-dump; do
  ln -sf "$(conda info --base)/envs/acomys-bio/bin/$t" \
         "$CONDA_PREFIX/bin/$t" 2>/dev/null
done

cp .env.example .env      # then add your NCBI API key (3x faster NCBI access)

python run_all.py --check        # preflight: env, tools, species gate
python run_all.py --dry-run      # show what each step would do
python run_all.py                # run everything
```

### Why two environments

The heavy bioinformatics tools pull their own Python, Perl and R pins. Solving
them alongside a pinned scientific stack is where conda environments go to die.
Split, both solve in under a minute.

### Setup trap worth knowing

**Never prepend `acomys-bio/bin` to `PATH`.** Conda envs contain full
interpreters, not just tools — prepending one hijacks `python`, and every
package check then silently describes the wrong environment. Use the symlinks
above, or append rather than prepend. `run_all.py --check` detects this
specific failure and names it.

Two related traps the preflight also catches:

- Conda does **not** error when a pinned package has no build for the resolved
  Python; it drops the package and reports success. A too-new Python therefore
  yields a half-populated env that looks fine.
- `conda activate` succeeding does not mean `python` comes from that env.
  Always confirm with `which python`.

After the first successful check, lock what actually resolved:

```bash
conda env export --no-builds > environment.lock.yml
```

---

## Code layout

```
code/
  config.py                        single source of truth - accessions, gene
                                   sets, seeds, thresholds, decision rules
  platform_compat.py               Windows/Linux tool resolution layer
  01_fetch_data.py                 data acquisition + provenance hashing
  02_build_ortholog_scaffold.py    orthology + region slicing   [LINCHPIN]
  02b_call_utrs.py                 UTR boundary calling
  03_kmer_phylo_controlled.py      alignment-free, phylogeny-controlled arm
  04_scrna_myeloid.py              single-cell myeloid composition
  05_selection_analysis.py         PAML branch-site models, multi-start
  06_expression_reanalysis.py      species-native re-quantification (kallisto)
  07_integration.py                sequence/expression convergence test
  08_figures.py                    all figures, generated from result CSVs
  09_test_nguyen_annotation.py     reading-frame integrity audit across three
                                   sources of Acomys coding sequence
  check_new_data.py                periodic re-check for newly public data
  inspect_alignments.py            alignment QC helper
  tests/
    make_synthetic_orthologs.py    positive/negative control fixtures
    debug_tier2.py                 region-recovery diagnostics
run_all.py                         single entry point
```

Numbered filenames are pipeline **steps**. `config.py` is not a step, hence no
number — a leading digit would make it un-importable.

---

## Design notes

**The k-mer arm is phylogenetically controlled.** *Acomys* is Deomyinae; *Mus*
is Murinae, and they split ~18.3 Ma. An Acomys-vs-Mus classifier reports
taxonomy, not phenotype, and will happily reach near-perfect accuracy meaning
nothing. The endpoint is therefore a **residual from phylogenetic expectation**
across ten species, compared against three matched control gene sets.
Deomyinae's sister clade is Gerbillinae, so *Meriones unguiculatus* — not
*Mus* — is the nearest annotated non-regenerative comparator.

**Expression is re-quantified species-natively.** The published kidney UUO
series quantified *Acomys* reads against the **mouse** transcriptome.
Cross-species pseudoalignment under-recovers divergent genes, which are
disproportionately chemokines and immune genes — exactly the targets here. The
published DE tables are used for orientation only.

**Sequence length is treated as a confound, not a nuisance.** It is estimated
within control genes only, where focal status cannot bias it, and both the
adjusted and unadjusted statistics are reported. Where they disagree, the
adjusted one is the answer. This rule is fixed in advance, in `config.py`.

---

## Carried design gaps

Run `python code/config.py` to print these. They live in the code so they
cannot be forgotten.

- **`gerbillinae_n1`** — Gerbillinae has one usable species. No clade-level
  claim about Gerbillinae is available.
- **`acomys_annotation_circularity`** — the *A. cahirinus* annotation is
  CAT-projected *from* GRCm39. Measuring divergence *from* mouse using
  mouse-projected gene models is circular. Analyses are therefore repeated
  under native, congener (*A. russatus*, independently RefSeq-annotated) and
  de novo transcript-assembly boundaries.
- **`clade_imbalance`** — Deomyinae 2, Gerbillinae 1, Murinae 5, Cricetidae 2.
  Leave-one-species-out sensitivity is required throughout.

---

## Working rules

1. **`config.py` is the only place constants live.** If a value appears twice,
   that is a bug.
2. **Raw data is never edited in place.** Downloads are hashed into
   `data/provenance_manifest.json`.
3. **No `_v2`, `_FIXED`, `_CORRECTED` files.** Use git.
4. **Figures are generated by script only.** No manual editing, and every
   plotted value is read from a result CSV rather than typed in.
5. **Decision rules are fixed before the analysis runs**, and are not relaxed
   after seeing results.
6. **Every accuracy claim carries a permutation null. Every effect size carries
   a bootstrap CI.**

---

## Verification

`code/tests/make_synthetic_orthologs.py` builds a synthetic ortholog scaffold
in exactly the layout the scaffold step writes and the k-mer step reads, with a
divergence signal planted in the focal module's 3'UTR and nowhere else. It ships
in two modes — a signal fixture the pipeline must detect and localise to the
correct region, and a null fixture in which it must find nothing anywhere.

```bash
cd code
python tests/make_synthetic_orthologs.py --clean --effect 0.30
python 03_kmer_phylo_controlled.py --arm residual --region utr3
python tests/make_synthetic_orthologs.py --clean --null      # must find nothing
```

Delete `data/reference/orthologs/` afterwards so synthetic and real data cannot
mix.

This is how one design bug was caught: an early composition null was
structurally incapable of firing, because shuffling every species independently
destroys homology and inflates all distances, so the observed value always sat
below the null. It was replaced with a composition-adjusted residual test, with
the shuffle retained only as a two-sided conservation diagnostic.

---

## Data

No sequencing reads, genome assemblies, transcriptomes or BLAST databases are
committed. All inputs are public and are re-downloaded by
`code/01_fetch_data.py`; accessions are declared in `code/config.py`.

Set `NCBI_API_KEY` in `.env` (see `.env.example`) before running the fetch step.

---

## Licence

MIT — see `LICENSE`.
