#!/usr/bin/env python3
"""
code_files.paths — Central configuration and sample discovery for the RNA-seq pipeline.

Design:
 - Read site values from environment variables (preferred).
 - Provide sensible defaults that mirror the previous hard-coded values.
 - Provide clear error messages when required folders/files are missing.
 - Expose the same names that existing scripts expect.

Usage:
  from code_files import paths
  print(paths.OUTDIR, len(paths.SAMPLES))
"""

from __future__ import annotations
import os
from pathlib import Path
from dataclasses import dataclass
from typing import List

# ----------------- small helper -----------------
def env_path(name: str, default: str | Path) -> Path:
    v = os.environ.get(name)
    if v:
        return Path(v)
    return Path(default)

# ----------------- Tool & root configuration (env-first) -----------------
PERLBASE = env_path("PERLBASE", "/sci/labs/zvika.granot/segev.munitz/softwares/perl5")
PSI_SIGMA = env_path("PSI_SIGMA", "/sci/labs/zvika.granot/segev.munitz/softwares/PSI-Sigma-2.1/dummyai.pl")
SALMON_BIN = env_path("SALMON_BIN", "/sci/labs/zvika.granot/segev.munitz/softwares/Salmon/salmon-latest_linux_x86_64/bin/salmon")

# Reference / index roots
REFERENCES_ROOT = env_path("REFERENCES_ROOT", "/sci/labs/zvika.granot/segev.munitz/references")
STAR_INDEX_ROOT = env_path("STAR_INDEX_ROOT", "/sci/labs/zvika.granot/segev.munitz/star_indexes")
SALMON_INDEX_ROOT = env_path("SALMON_INDEX_ROOT", "/sci/labs/zvika.granot/segev.munitz/salmon_indexes")

# FASTQ suffixes (defaults)
R1_SUFFIX = os.environ.get("R1_SUFFIX", "_1.fastq.gz")
R2_SUFFIX = os.environ.get("R2_SUFFIX", "_2.fastq.gz")

@dataclass(frozen=True)
class Comparison:
    name: str
    groupa: Path
    groupb: Path

# ----------------- Run-specific (edit via env or configs/site.env) -----------------
NAME = os.environ.get("RUN_NAME", os.environ.get("NAME", "run_1"))
OUTDIR = Path(os.environ.get("OUTDIR", "/sci/labs/zvika.granot/segev.munitz/psi_sigma_outputs")) / NAME
FASTQ_DIR = Path(os.environ.get("FASTQ_DIR", str(OUTDIR / "gz_files")))

KIND = os.environ.get("KIND", "Mouse")        # "Human" or "Mouse"
VERSION = os.environ.get("VERSION", "M38")    # e.g. "M38" or "v45"
READ_LENGTH = int(os.environ.get("READ_LENGTH", "100"))
sjdbOverhang = max(READ_LENGTH - 1, 1)

# Derived genome prefix (used in FASTA filename patterns)
GENOME_PRE_FIX = "GRCh38" if KIND == "Human" else "GRCm39"

# ----------------- GENCODE/URLs (defaults, can be overridden) -----------------
if KIND == "Human":
    GENCODE_BASE_URL = os.environ.get("GENCODE_BASE_URL", "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_45")
else:
    GENCODE_BASE_URL = os.environ.get("GENCODE_BASE_URL", "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M38")

# Provide optional download URLs (used by helper scripts if needed)
GENOME_URL = os.environ.get("GENOME_URL", "")
ANNOTATION_URL = os.environ.get("ANNOTATION_URL", "")
TRANSCRIPTS_URL = os.environ.get("TRANSCRIPTS_URL", "")

# ----------------- Reference locations derived from roots -----------------
REF_BUNDLE_DIR = Path(REFERENCES_ROOT) / KIND / f"Genecode_{VERSION}_rl{READ_LENGTH}"
STAR_INDEX_DIR = Path(STAR_INDEX_ROOT) / KIND / f"Genecode_{VERSION}_rl{READ_LENGTH}"
SALMON_INDEX_DIR = Path(SALMON_INDEX_ROOT) / KIND / f"Genecode_{VERSION}_rl{READ_LENGTH}"

GENOME_FASTA = REF_BUNDLE_DIR / f"{GENOME_PRE_FIX}.primary_assembly.genome.fa"
ANNOTATION_GTF = REF_BUNDLE_DIR / f"gencode.{VERSION}.annotation.gtf"
SALMON_TRANSCRIPTS_FASTA = REF_BUNDLE_DIR / f"gencode.{VERSION}.transcripts.fa"

# ----------------- Sample discovery -----------------
def build_samples_from_fastq_dir(
    fastq_dir: Path,
    r1_suffix: str = R1_SUFFIX,
    r2_suffix: str = R2_SUFFIX,
    require_paired: bool = True,
) -> List[dict]:
    fastq_dir = Path(fastq_dir)
    if not fastq_dir.exists() or not fastq_dir.is_dir():
        raise SystemExit(f"ERROR: FASTQ_DIR does not exist or is not a directory: {fastq_dir}\n"
                         f"Set env var FASTQ_DIR or update paths.py. Current OUTDIR={OUTDIR}")
    samples = []
    r1_files = sorted(fastq_dir.glob(f"*{r1_suffix}"))
    if not r1_files:
        raise SystemExit(f"ERROR: No R1 files found in {fastq_dir} matching pattern '*{r1_suffix}'.")
    for r1 in r1_files:
        sample_id = r1.name[:-len(r1_suffix)]
        r2 = fastq_dir / f"{sample_id}{r2_suffix}"
        if require_paired:
            if not r2.exists():
                # skip incomplete pairs but warn
                continue
            samples.append({"id": sample_id, "r1": r1.resolve(), "r2": r2.resolve()})
        else:
            entry = {"id": sample_id, "r1": r1.resolve()}
            if r2.exists():
                entry["r2"] = r2.resolve()
            samples.append(entry)
    if require_paired and not samples:
        raise SystemExit(f"ERROR: Found R1 files but no complete R1/R2 pairs using "
                         f"r1='{r1_suffix}', r2='{r2_suffix}' in {fastq_dir}.")
    samples.sort(key=lambda d: str(d["id"]))  # type: ignore[arg-type]
    return samples

# Build SAMPLES now (will exit early if FASTQ_DIR is missing or empty)
SAMPLES = build_samples_from_fastq_dir(FASTQ_DIR)

# ----------------- Computation defaults and thresholds -----------------
THREADS = int(os.environ.get("THREADS", "16"))
DO_QC = os.environ.get("DO_QC", "False").lower() in ("1", "true", "yes")
DO_TRIM = os.environ.get("DO_TRIM", "False").lower() in ("1", "true", "yes")

SALMON_LIBTYPE = os.environ.get("SALMON_LIBTYPE", "A")
SALMON_EXTRA_ARGS = os.environ.get("SALMON_EXTRA_ARGS", "")
SALMON_TPM_THRESHOLD = float(os.environ.get("SALMON_TPM_THRESHOLD", "5.0"))
SALMON_FILTER_MODE = os.environ.get("SALMON_FILTER_MODE", "either")  # "either" or "both"

PSISIGMA_ABSPSI_MIN = float(os.environ.get("PSISIGMA_ABSPSI_MIN", "20.0"))
PSISIGMA_P_MAX = float(os.environ.get("PSISIGMA_P_MAX", "0.05"))
PSISIGMA_FDR_MAX = float(os.environ.get("PSISIGMA_FDR_MAX", "0.05"))

# Default group file paths (these will be created by step3)
GROUP_HEALTHY = OUTDIR / "groups" / "H.bams.txt"
GROUP_SICK_1 = OUTDIR / "groups" / "N1.bams.txt"
GROUP_SICK_2 = OUTDIR / "groups" / "N2.bams.txt"

COMPARISONS = [
    Comparison("H_vs_N1", GROUP_HEALTHY, GROUP_SICK_1),
    Comparison("H_vs_N2", GROUP_HEALTHY, GROUP_SICK_2),
    Comparison("N1_vs_N2", GROUP_SICK_1, GROUP_SICK_2),
]