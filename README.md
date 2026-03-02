# RNA-seq Pipeline for SLURM and Nextflow

Comprehensive pipeline for bulk RNA-seq processing with:

- **STAR** alignment
- **featureCounts** quantification
- **PSI-Sigma** differential splicing
- **Salmon** transcript quantification
- **Salmon-based filtering** of PSI-Sigma outputs

This repository supports **two orchestration modes**:

1. **Nextflow-first (recommended):** `main.nf` + `nextflow.config`
2. **Native SLURM mode:** `sbatch-commands/*.slurm` + `submit_all.sh`

The Python scripts in `code-files/` are the core implementation and remain the source of truth for workflow behavior.

---

## Table of Contents
1. [Pipeline Overview](#pipeline-overview)
2. [Execution Modes](#execution-modes)
3. [Prerequisites & Environment](#prerequisites--environment)
4. [Configuration in `paths.py`](#configuration-in-pathspy)
5. [Input Requirements](#input-requirements)
6. [Run with Nextflow (recommended)](#run-with-nextflow-recommended)
7. [Run with SLURM scripts](#run-with-slurm-scripts)
8. [Step-by-Step Detailed Behavior](#step-by-step-detailed-behavior)
9. [Dependency Graph](#dependency-graph)
10. [Output Directory Structure](#output-directory-structure)
11. [Reference/Index Directory Structure](#referenceindex-directory-structure)
12. [Troubleshooting](#troubleshooting)
13. [Operational Notes & Best Practices](#operational-notes--best-practices)

---

## Pipeline Overview

The pipeline is organized into seven stages:

| Step | Script | What it does |
|---|---|---|
| 1 | `step1_prep_run.py` | Initializes run directory, writes `samples.txt`, builds/reuses STAR and Salmon indexes |
| 2 | `step2_star_align_core.py` | STAR alignment per sample (array-style), produces sorted BAM + BAI |
| 3 | `step3_make_group_files.py` | Creates group BAM list files from sample/group mapping |
| 4 | `step4_featurecounts.py` | Runs featureCounts and creates cleaned gene-count matrix |
| 5 | `step5_run_psi_sigma.py` | Runs PSI-Sigma comparisons from `COMPARISONS` |
| 6 | `step6_salmon_quant_core.py` | Salmon quantification per sample (array-style) |
| 7 | `step7_psi_sigma_filtered_by_salmon.py` | Filters PSI-Sigma events by Salmon TPM and significance thresholds |

---

## Execution Modes

### 1) Nextflow mode (recommended)

Use when you want:
- resumability (`-resume`)
- cleaner DAG-based orchestration
- easier resource management in one place
- less manual dependency chaining

Key files:
- `main.nf`
- `nextflow.config`

### 2) SLURM script mode

Use when you want:
- direct/manual `sbatch` control
- one script per stage
- compatibility with existing team habits

Key files:
- `sbatch-commands/sbatch_01_prep.slurm` ... `sbatch_07_salmon_filter.slurm`
- `sbatch-commands/submit_all.sh`

---

## Prerequisites & Environment

Ensure these tools are available on compute nodes (in `PATH` or via module/conda):

### Required software
- Python 3.8+
- STAR (2.7.x recommended)
- Samtools
- Subread (`featureCounts`)
- Salmon
- Perl (with PSI-Sigma dependencies, e.g. `PDL`, `Statistics::R`)
- PSI-Sigma (`dummyai.pl`)
- Nextflow (for Nextflow mode)

### Python packages
```bash
pip install pandas openpyxl --break-system-packages
```

### Optional environment checks
```bash
python3 --version
STAR --version
samtools --version
featureCounts -v
salmon --version
perl -v
nextflow -version
```

---

## Configuration in `paths.py`

Single source of configuration:

- `code-files/paths.py`

### Primary fields to edit
- **Run naming/paths**
  - `NAME`
  - `OUTDIR`
  - `FASTQ_DIR`
- **Reference selection**
  - `KIND` (`Human` / `Mouse`)
  - `VERSION`
  - `READ_LENGTH`
- **Tool locations**
  - `PERLBASE`
  - `PSI_SIGMA`
  - `SALMON_BIN`
- **Roots for references/indexes**
  - `REFERENCES_ROOT`
  - `STAR_INDEX_ROOT`
  - `SALMON_INDEX_ROOT`
- **PSI-Sigma comparisons**
  - `COMPARISONS`
  - related group files
- **Filtering thresholds**
  - `SALMON_TPM_THRESHOLD`
  - `SALMON_FILTER_MODE`
  - `PSISIGMA_ABSPSI_MIN`, `PSISIGMA_P_MAX`, `PSISIGMA_FDR_MAX`

### FASTQ naming conventions
Defaults:
```python
R1_SUFFIX = "_1.fastq.gz"
R2_SUFFIX = "_2.fastq.gz"
```

### Example (trimmed)
```python
NAME = "run_1"
OUTDIR = Path("/path/to/output") / NAME
FASTQ_DIR = OUTDIR / "gz_files"

KIND = "Mouse"
VERSION = "M38"
READ_LENGTH = 100

GROUP_HEALTHY = OUTDIR / "groups" / "H.bams.txt"
GROUP_SICK_1  = OUTDIR / "groups" / "N1.bams.txt"
GROUP_SICK_2  = OUTDIR / "groups" / "N2.bams.txt"

COMPARISONS = [
    Comparison("H_vs_N1", GROUP_HEALTHY, GROUP_SICK_1),
    Comparison("H_vs_N2", GROUP_HEALTHY, GROUP_SICK_2),
    Comparison("N1_vs_N2", GROUP_SICK_1, GROUP_SICK_2),
]
```

> **Important:** Keep `COMPARISONS` group file paths consistent with where Step 3 writes group files in your configuration/layout.

---

## Input Requirements

### 1) FASTQ files
- Must exist in `FASTQ_DIR`
- Must match suffixes (`R1_SUFFIX`, `R2_SUFFIX`)
- For paired-end mode, complete R1/R2 pairs are required

### 2) Sample-to-group mapping
Create a `groups.tsv` before Step 3.

Current Step 3 code expects:
- `<OUTDIR>/star_alignments/groups.tsv`

Format (no header):
```tsv
sample1	H
sample2	H
sample3	N1
sample4	N1
sample5	N2
sample6	N2
```

### 3) Reference bundles/indexes
- STAR and Salmon indexes are built/reused automatically in Step 1
- Missing FASTA/GTF/transcript files can be downloaded (as configured)

---

## Run with Nextflow (recommended)

### Files
- `main.nf`: full DAG orchestration
- `nextflow.config`: resources, queue, concurrency, profiles

### Default process resources (from `nextflow.config`)
- Queue: `glacier`
- `PREP`: 12 CPUs, 60 GB, 24h
- `STAR_ALIGN`: 12 CPUs, 60 GB, 24h, `maxForks = 8`
- `MAKE_GROUPS`: 1 CPU, 2 GB, 10m
- `FEATURE_COUNTS`: 12 CPUs, 60 GB, 24h
- `PSI_SIGMA`: 16 CPUs, 60 GB, 24h
- `SALMON_QUANT`: 16 CPUs, 32 GB, 24h, `maxForks = 8`
- `SALMON_FILTER`: 4 CPUs, 16 GB, 24h

### Basic commands
```bash
# From repo root
nextflow run main.nf -profile slurm

# Resume interrupted run
nextflow run main.nf -profile slurm -resume
```

### Useful parameter overrides
```bash
# Use a different python executable
nextflow run main.nf -profile slurm --python python3.11

# Disable featureCounts paired-end mode
nextflow run main.nf -profile slurm --featurecounts_paired_end false
```

### Profiles
- `-profile slurm`: use SLURM executor
- `-profile standard`: local execution (useful for debugging small tests)

### What Nextflow is doing internally
- Runs Step 1 prep and reads `<OUTDIR>/samples.txt`
- Builds per-sample task tuples
- Fans out STAR and Salmon per sample
- Preserves stage ordering equivalent to `submit_all.sh`
- Exports `SLURM_ARRAY_TASK_ID` per sample for script compatibility

### Nextflow runtime artifacts
- `work/` task working directories
- `.nextflow/` metadata
- `.nextflow.log` run log

---

## Run with SLURM scripts

### Option A: run full chain automatically
```bash
cd sbatch-commands
chmod +x submit_all.sh
./submit_all.sh
```

### Option B: submit manually with dependencies
```bash
# 1) Prep
jid_prep=$(sbatch sbatch_01_prep.slurm | awk '{print $4}')

# 2) Count samples
N=$(grep -v '^\s*$' samples.txt | wc -l)

# 3) STAR after prep
jid_star=$(sbatch --dependency=afterok:${jid_prep} --array=1-${N}%8 sbatch_02_star_array.slurm | awk '{print $4}')

# 4) Groups + featureCounts after STAR
jid_groups=$(sbatch --dependency=afterok:${jid_star} sbatch_03_make_groups.slurm | awk '{print $4}')
jid_fc=$(sbatch --dependency=afterok:${jid_star} sbatch_04_featurecounts.slurm | awk '{print $4}')

# 5) PSI-Sigma after groups + featureCounts
jid_psisigma=$(sbatch --dependency=afterok:${jid_groups}:${jid_fc} sbatch_05_run_psi_sigma.slurm | awk '{print $4}')

# 6) Salmon after STAR
jid_salmon=$(sbatch --dependency=afterok:${jid_star} --array=1-${N}%8 sbatch_06_salmon_array.slurm | awk '{print $4}')

# 7) Final filter after PSI-Sigma + Salmon
sbatch --dependency=afterok:${jid_psisigma}:${jid_salmon} sbatch_07_salmon_filter.slurm
```

---

## Step-by-Step Detailed Behavior

### Step 1: Prep (`step1_prep_run.py`)
- Creates `OUTDIR`
- Writes `samples.txt` from discovered FASTQ pairs
- Builds/reuses STAR index
- Builds/reuses Salmon index
- Optional reference downloads when missing

### Step 2: STAR (`step2_star_align_core.py`)
- One task per sample
- Resolves sample by task ID (`--task-id` or `SLURM_ARRAY_TASK_ID`)
- Validates FASTQ files
- Runs STAR and writes sorted BAM
- Runs `samtools index`

### Step 3: Groups (`step3_make_group_files.py`)
- Reads sample IDs from `samples.txt`
- Reads mapping from `groups.tsv`
- Validates BAM paths for mapped samples
- Writes `<group>.bams.txt` files for PSI-Sigma input

### Step 4: featureCounts (`step4_featurecounts.py`)
- Uses explicit `--bam` inputs or auto-discovers STAR BAMs
- Runs `featureCounts`
- Produces:
  - `featureCounts.gene_counts.txt` (raw)
  - `gene_counts.matrix.tsv` (clean)

### Step 5: PSI-Sigma (`step5_run_psi_sigma.py`)
- Validates GTF, PSI-Sigma script, and group files
- Sets Perl environment (`PATH`, `PERL5LIB`)
- Runs one PSI-Sigma job per comparison in `COMPARISONS`

### Step 6: Salmon (`step6_salmon_quant_core.py`)
- One task per sample
- Ensures Salmon index exists
- Task 1 can build missing index; others wait
- Runs Salmon quant and writes `quant.sf`

### Step 7: Salmon filter (`step7_psi_sigma_filtered_by_salmon.py`)
- Loads PSI-Sigma outputs per comparison
- Loads group-specific Salmon TPMs
- Applies TPM + statistical filters
- Writes filtered `.tsv` and `.xlsx`

---

## Dependency Graph

```text
Prep (1)
  ↓
STAR Array (2)
  ↓
  ├─→ Groups (3) ──→ PSI-Sigma (5) ──┐
  ├─→ featureCounts (4) ─────────────┤
  └─→ Salmon Array (6) ──────────────┤
                                      ↓
                              Filter Results (7)
```

---

## Output Directory Structure

```text
<OUTDIR>/
├── gz_files/                                  # FASTQ input directory (default FASTQ_DIR)
│   ├── <sample>_1.fastq.gz
│   └── <sample>_2.fastq.gz
│
├── samples.txt                                # list of sample IDs (Step 1)
│
├── star_alignments/
│   ├── groups.tsv                             # sample->group mapping (user-provided for Step 3)
│   ├── bam/                                   # STAR outputs (Step 2)
│   │   └── <sample>/
│   │       ├── <sample>.Aligned.sortedByCoord.out.bam
│   │       ├── <sample>.Aligned.sortedByCoord.out.bam.bai
│   │       ├── <sample>.Log.final.out
│   │       ├── <sample>.Log.out
│   │       ├── <sample>.Log.progress.out
│   │       └── <sample>.SJ.out.tab
│   ├── logs/
│   ├── tmp/
│   └── groups/                                # group BAM lists (Step 3)
│       ├── H.bams.txt
│       ├── N1.bams.txt
│       └── N2.bams.txt
│
├── counts/                                    # Step 4
│   ├── featureCounts.gene_counts.txt
│   └── gene_counts.matrix.tsv
│
├── salmon/                                    # Step 6
│   └── <sample>/
│       ├── quant.sf
│       ├── quant.genes.sf
│       ├── cmd_info.json
│       ├── lib_format_counts.json
│       ├── logs/
│       └── aux_info/
│
├── <comparison_name>/                         # Step 5 + Step 7 outputs
│   ├── *.PSIsigma*.txt                        # unfiltered PSI-Sigma
│   ├── *.salmon_tpm*.tsv                      # filtered
│   └── *.salmon_tpm*.xlsx                     # filtered
│
├── logs/
└── tmp/
```

---

## Reference/Index Directory Structure

Configured in `paths.py`:

```text
<REFERENCES_ROOT>/<KIND>/Genecode_<VERSION>_rl<READ_LENGTH>/
├── <GENOME_PRE_FIX>.primary_assembly.genome.fa
├── gencode.<VERSION>.annotation.gtf
└── gencode.<VERSION>.transcripts.fa

<STAR_INDEX_ROOT>/<KIND>/Genecode_<VERSION>_rl<READ_LENGTH>/
├── Genome
├── SA
├── SAindex
└── ...

<SALMON_INDEX_ROOT>/<KIND>/Genecode_<VERSION>_rl<READ_LENGTH>/
├── info.json
├── ctable.bin
└── ...
```

---

## Troubleshooting

### Nextflow not available
```bash
nextflow -version
```
If missing, install Nextflow and re-run.

### No samples found
- Check `FASTQ_DIR`.
- Verify `R1_SUFFIX`/`R2_SUFFIX` match filenames exactly.
- Ensure paired files exist when `require_paired=True`.

### STAR index build failures
- Increase memory/time for prep step.
- Verify genome FASTA and GTF paths/URLs.
- Confirm files are valid and readable.

### Step 2 failures (STAR)
- Confirm STAR and samtools are in `PATH`.
- Check `SLURM_ARRAY_TASK_ID` / `samples.txt` line mapping.
- Validate FASTQ paths in discovered samples.

### Step 3 fails to find `groups.tsv`
- Ensure mapping file exists at the configured/expected location.
- Ensure sample IDs match `samples.txt` IDs.

### featureCounts reports no BAMs
- Auto-discovery path: `OUTDIR/star_alignments/bam/*/*.Aligned.sortedByCoord.out.bam`
- Confirm Step 2 completed and path is correct.

### PSI-Sigma Perl/module errors
- Verify `PERLBASE` and `PSI_SIGMA` in `paths.py`.
- Check required Perl modules (e.g., `PDL`, `Statistics::R`).

### Salmon quant/index issues
- Validate `SALMON_BIN` path.
- Verify index directory content (`info.json`, `ctable.bin`).
- Waiting behavior for index build in arrays is expected.

### Empty Step 7 filtered results
- Consider lowering `SALMON_TPM_THRESHOLD`.
- Try `SALMON_FILTER_MODE = "either"`.
- Confirm PSI-Sigma transcript IDs and Salmon transcript IDs match sufficiently.

### Resource/time failures
- Increase CPU/memory/time in SLURM scripts or `nextflow.config`.
- Reduce concurrency (`%N` in arrays / `maxForks` in Nextflow).

### Permission errors
- Ensure write access to `OUTDIR`, reference/index roots, and log dirs.
- Check quotas and free disk (`df -h`, `quota -s`).

---

## Operational Notes & Best Practices

- Prefer **Nextflow + `-resume`** for robust long-running orchestration.
- Keep all run-specific settings in `paths.py`; avoid hardcoding in scripts.
- Pin tool versions in your environment/modules for reproducibility.
- Keep a dedicated per-run `OUTDIR` (`NAME`) to avoid accidental overwrite.
- Validate `groups.tsv` and `COMPARISONS` before expensive compute steps.
- Inspect logs after each stage before launching downstream analysis.

---

**Pipeline Version:** 2.1 
**Last Updated:** February 2025

