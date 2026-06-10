# =============================================================================
# MutSeqR BaP Reproduction — Project 1
# Script 06: COSMIC Signature Fitting
# =============================================================================
#
# OVERVIEW
# --------
# This script reproduces Section 3.5 (Signature Assignment) of:
#   Dodge et al. (2025) "MutSeqR: an open source R package for standardized
#   analysis of error-corrected next-generation sequencing data in genetic
#   toxicology". Bioinformatics Advances, vbaf265.
#
# signature_fitting() uses SigProfilerAssignment (Python) to assign COSMIC
# SBS signatures to the 96-base trinucleotide mutation spectrum of each
# dose group. Mutations are aggregated by group and compared against the
# COSMIC signature database to identify which known mutational processes
# contribute to the observed spectrum.
#
# BaP's expected signature is SBS4 — associated with tobacco smoke exposure
# and bulky DNA adducts. BaP is a major component of tobacco smoke, so a
# dominant SBS4 contribution confirms the assay is detecting BaP's known
# mechanism.
#
# PYTHON ENVIRONMENT NOTE
# ------------------------
# signature_fitting() creates a virtual environment via reticulate. In this
# session, reticulate auto-installed SigProfilerAssignment and
# SigProfilerMatrixGenerator into an ephemeral environment (Python 3.12).
# We use python_version = "3.12" to match.
#
# OUTPUT STRUCTURE
# -----------------
# Results are written to disk (not returned as R objects). The key output is:
# outputs/SigProfiler/{project_name}/matrices/output/Assignment_Solution/
#   Activities/ — signature contribution percentages per group
#   SampleReconstruction/WebPNGs/ — visual summaries per group
#
# Cosine similarity > 0.9 indicates a robust reconstruction.
#
# REPRODUCTION TARGETS (section 3.5)
# ------------------------------------
# - SBS4 dominant at all BaP dose groups (dose-dependent increase)
# - SBS4 contribution at HIGH dose: 97.51%
# - Additional minor signatures: SBS1, SBS5 (aging), SBS94
# - Cosine similarity > 0.9 for all dose groups
#
# =============================================================================

library(MutSeqR)
library(ExperimentHub)
library(dplyr)

cat("=== Script 06: COSMIC Signature Fitting ===\n\n")

# ---- 1. Load filtered mutation data -----------------------------------------
# signature_fitting() requires the raw filtered mutation data (not MF
# summaries). It filters to SNVs internally and uses the filter_mut column
# to exclude flagged variants.

bap_filtered <- readRDS("outputs/01_bap_filtered.rds")
cat("Loaded filtered data:", nrow(bap_filtered), "rows\n")
cat("filter_mut column present:",
    "filter_mut" %in% colnames(bap_filtered), "\n\n")

# ---- 2. Set up output directory ---------------------------------------------

dir.create("outputs/SigProfiler", showWarnings = FALSE, recursive = TRUE)
output_path <- file.path(getwd(), "outputs")
cat("Output path:", output_path, "\n\n")

# ---- 3. Run signature fitting -----------------------------------------------
# group = "dose_group": aggregate mutations by dose group (Control, Low,
#   Medium, High) — matching the paper's approach of grouping by dose.
# project_genome = "mm10": MutaMouse uses the mm10 reference genome.
# python_version = "3.12": matches the reticulate auto-installed version.
#
# NOTE: On first run, this will download the mm10 reference genome for
# SigProfilerMatrixGenerator (~500 MB). This is a one-time download.
# Subsequent runs will use the cached genome.
#
# Runtime: several minutes on first run (genome download + signature fitting).

cat("=== Running signature_fitting() ===\n")
cat("Grouping by: dose_group\n")
cat("Reference genome: mm10\n")
cat("Python version: 3.12\n")
cat("NOTE: First run downloads mm10 genome (~500 MB)\n\n")

tryCatch({
  signature_fitting(
    mutation_data  = bap_filtered,
    project_name   = "BaP_reproduction",
    project_genome = "mm10",
    env_name       = "MutSeqR",
    group          = "dose_group",
    output_path    = output_path,
    python_version = "3.12"
  )
  cat("\nsignature_fitting() completed successfully\n")
}, error = function(e) {
  cat("\nERROR in signature_fitting():\n")
  cat(conditionMessage(e), "\n")
  cat("\nThis may be due to Python environment conflicts.\n")
  cat("Try running setup_mutseqr_python() first:\n")
  cat("  setup_mutseqr_python(force = TRUE)\n")
})

# ---- 4. Read and display results --------------------------------------------
# Results are written to disk. We read the Activities file to get
# signature contribution percentages per dose group.

cat("\n=== Reading signature assignment results ===\n")

# The activities file path follows SigProfilerAssignment's output structure
activities_path <- file.path(
  output_path, "SigProfiler", "BaP_reproduction",
  "matrices", "output", "Assignment_Solution",
  "Activities", "SBS96_Assignment_Solution_Activities.txt"
)

if (file.exists(activities_path)) {
  activities <- read.delim(activities_path, check.names = FALSE)
  cat("Activities file found\n")
  cat("Dimensions:", nrow(activities), "rows x", ncol(activities), "cols\n\n")
  print(activities)
} else {
  cat("Activities file not found at expected path.\n")
  cat("Expected:", activities_path, "\n")
  cat("Checking outputs/SigProfiler for actual structure:\n")
  
  if (dir.exists(file.path(output_path, "SigProfiler"))) {
    # List all files to find the actual output
    all_files <- list.files(
      file.path(output_path, "SigProfiler"),
      recursive = TRUE, full.names = FALSE
    )
    cat("\nFiles in SigProfiler output:\n")
    print(head(all_files, 30))
  } else {
    cat("SigProfiler output directory not created yet.\n")
    cat("Check for errors above.\n")
  }
}

# ---- 5. Reproduction check --------------------------------------------------

cat("\n=== REPRODUCTION CHECK ===\n")
cat("Target: SBS4 dominant at all BaP dose groups\n")
cat("Target: SBS4 contribution at High dose = 97.51%\n")
cat("Target: Cosine similarity > 0.9 for all groups\n\n")

# Check if we can find the SBS4 contributions
if (exists("activities") && !is.null(activities)) {
  if ("SBS4" %in% colnames(activities)) {
    cat("SBS4 contributions by dose group:\n")
    sbs4_cols <- activities[, c(1, which(colnames(activities) == "SBS4"))]
    print(sbs4_cols)
  } else {
    cat("SBS4 column not found in activities. Columns present:\n")
    print(colnames(activities))
  }
}

cat("\nNote: Full visual results (reconstructed profiles, signature\n")
cat("contributions) are in:\n")
cat(file.path(output_path, "SigProfiler", "BaP_reproduction",
              "matrices", "output", "Assignment_Solution",
              "SampleReconstruction", "WebPNGs"), "\n")

# =============================================================================
# END OF SCRIPT 06
# =============================================================================