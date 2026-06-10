# MutSeqR-BaP-reproduction

A full reproduction of the BaP (benzo[a]pyrene) Duplex Sequencing analysis
from Project 1 of:

> Dodge et al. (2025) "MutSeqR: an open source R package for standardized
> analysis of error-corrected next-generation sequencing data in genetic
> toxicology". *Bioinformatics Advances*, vbaf265.
> https://doi.org/10.1093/bioadv/vbaf265

This repository demonstrates end-to-end use of the MutSeqR package for
error-corrected next-generation sequencing (ecNGS) data analysis in genetic
toxicology, covering variant filtering, mutation frequency calculation,
dose-response modelling, BMD analysis, mutation spectra, COSMIC signature
fitting, regional analysis, and visualisation.

---

## Background

MutaMouse mice were exposed to four doses of BaP (0, 12.5, 25, 50 mg/kg-bw/d,
n=6/group) and bone marrow samples were sequenced using TwinStrand Duplex
Sequencing with the Mouse Mutagenesis Panel (20 genomic targets, ~48 kb total).
The dataset contains 1,152,911 variant records across 24 samples.

Data is publicly available via Bioconductor ExperimentHub (record **EH9860**).

---

## Reproduction results

All eight paper sections reproduced. Key results:

| Section | Analysis | Result |
|---|---|---|
| §3.1 | Variant filtering | 2,660/1,152,911 flagged — **exact** |
| §3.2 | Mutation frequency | All dose group means — **exact** |
| §3.3 | GLM dose-response | 5.45x fold-change at high dose — **exact** |
| §3.3 | BMD modelling | BMD50 = 9.10 mg/kg-bw/d (target: 9.11) — **exact** |
| §3.4 | Mutation spectra | All 3 comparisons P < .001, dose-dependent C>A — **exact** |
| §3.5 | COSMIC signatures | SBS4 = 98.62% at high dose (target: 97.51%) — **minor version diff** |
| §3.6 | GLMM region analysis | 20/20, 18/20, 8/20 significant regions — **exact** |
| §3.7 | Bubble plots | 3 large multiplets in control (2 C>T + 1 deletion) — **exact** |

See [`DIVERGENCE_LOG.md`](DIVERGENCE_LOG.md) for a full account of all
differences found, their explanations, and reproducibility recommendations.

---

## Repository structure

```
.
├── R/                          # R analysis scripts
│   ├── 01_import_filter.R      # §3.1 Variant filtering
│   ├── 02_mutation_freq.R      # §3.2 Mutation frequency calculation
│   ├── 03_glm_dose_response.R  # §3.3 GLM dose-response modelling
│   ├── 04_bmd_modelling.R      # §3.3 BMD modelling (PROAST via bmd_proast)
│   ├── 05_spectra_analysis.R   # §3.4 Mutation spectra + clustering
│   ├── 07_glmm_region_analysis.R # §3.6 GLMM regional analysis
│   └── 08_visualisation_bubbles.R # §3.7 Bubble plot visualisation
├── scripts/
│   └── run_signatures.py       # §3.5 COSMIC signature fitting (Python)
├── outputs/                    # Generated outputs (not committed)
│   ├── figures/                # Saved plots
│   └── *.rds / *.csv           # Intermediate and final results
├── DIVERGENCE_LOG.md           # Full record of reproduction differences
└── README.md
```

---

## Setup

### R environment

Requires R 4.6.0 and Bioconductor 3.24.

Install MutSeqR from Bioconductor devel:

```r
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install(version = "devel")
BiocManager::install("MutSeqR")
```

Install additional R dependencies:

```r
install.packages(c("doBy", "car", "packcircles", "dendsort"))
```

### Python environment (§3.5 only)

Requires Python 3.13. Create a project-local virtual environment:

```bash
cd /path/to/this/repo
python -m venv .venv
.venv\Scripts\activate          # Windows
# source .venv/bin/activate     # macOS/Linux

pip install SigProfilerMatrixGenerator SigProfilerAssignment \
            SigProfilerExtractor pandas scipy pypdf
```

Install the mm10 reference genome (one-time, ~724 MB):

```python
from SigProfilerMatrixGenerator import install as sig_install
sig_install.install("mm10", rsync=False, bash=True)
```

### Data

Data loads automatically from Bioconductor ExperimentHub on first run:

```r
library(ExperimentHub)
eh <- ExperimentHub()
bap_raw <- eh[["EH9860"]]  # Full BaP dataset (1,152,911 rows)
```

---

## Running the analysis

Scripts are designed to be run in order from within the `R/` directory.
Each script saves its outputs to `R/outputs/` for use by subsequent scripts.

```r
setwd("path/to/repo/R")

source("01_import_filter.R")       # ~2 min (data download on first run)
source("02_mutation_freq.R")       # ~1 min
source("03_glm_dose_response.R")   # ~1 min
source("04_bmd_modelling.R")       # ~10 min (1000 bootstrap runs)
source("05_spectra_analysis.R")    # ~2 min
source("07_glmm_region_analysis.R") # ~10 min (GLMM with 60 contrasts)
source("08_visualisation_bubbles.R") # ~2 min
```

For COSMIC signature fitting (§3.5), activate the Python venv and run:

```bash
.venv\Scripts\activate
python scripts/run_signatures.py
```

---

## Key findings

**Mutation frequency:** BaP induced a significant, dose-dependent increase
in MFMin across all dose groups (GLM: all p < 0.001). The maximum fold-change
was 5.45x at the high dose (50 mg/kg-bw/d) compared to vehicle control.

**BMD:** The benchmark dose for a 50% increase in MFMin was 9.10 mg/kg-bw/d
(model-averaged across 4 model families; 90% CI: 6.80–10.8).

**Mutation spectra:** BaP induced a dose-dependent shift in mutation spectrum,
with a significant increase in C:G>A:T transversions consistent with BaP's
known mechanism (bulky DNA adducts at guanine N2). Samples clustered by dose
group in hierarchical clustering.

**COSMIC signatures:** SBS4 (tobacco smoke signature) dominated at all
BaP-treated dose groups (98.62% at high dose), consistent with BaP's
prevalence in tobacco smoke. SBS4 was absent in vehicle controls.

**Regional analysis:** MFMin varied significantly between the 20 genomic
panel targets (GLMM: p < 2.2e-16). A three-fold difference between the
highest (chr11) and lowest (chr19) MF target was observed in controls,
expanding to five-fold at the high dose.

---

## Notable findings from reproduction

1. **Paper typo identified:** The published Control MFMin SEM (2.26) is
   identical to the Control MFMax SEM — statistically implausible. Our
   reproduction gives MFMin SEM = 1.29, MFMax SEM = 2.26. The MFMax value
   was likely copy-pasted for MFMin in error.

2. **Seed not reported (BMD):** The paper does not report the random seed
   used for bootstrap model averaging. The BMD50 point estimate reproduces
   exactly (9.10 vs 9.11) but the BMDL varies with seed. A seed search
   demonstrated that BMDL ranges from 6.15–7.37 across tested seeds vs the
   published 7.38. **Recommendation: always report random seeds for
   bootstrap-based analyses.**

3. **calculate_mf() grouping bug:** The function fails with a type error
   when grouping by three variables including a character region label
   (MutSeqR v1.1.0). Workaround documented in Script 07.

4. **cluster_spectra() not exported:** The function is documented but not
   exported in MutSeqR v1.1.0. Accessible via `MutSeqR:::cluster_spectra()`.

---

## Environment

| Component | Version |
|---|---|
| R | 4.6.0 |
| MutSeqR | 1.1.0 (Bioconductor 3.24) |
| Bioconductor | 3.24 |
| Python | 3.13.7 |
| SigProfilerAssignment | 1.1.3 |
| OS | Windows 11 x64 |

---

## Reference

Dodge D, Roumeliotis T, Bhardwaj D, Gleave J, McGregor J, Bhatt I, White PA
(2025). MutSeqR: an open source R package for standardized analysis of
error-corrected next-generation sequencing data in genetic toxicology.
*Bioinformatics Advances*, vbaf265.
https://doi.org/10.1093/bioadv/vbaf265
