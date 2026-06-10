# =============================================================================
# MutSeqR BaP Reproduction — Project 1
# Script 05: Mutation Spectra Analysis
# =============================================================================
#
# OVERVIEW
# --------
# This script reproduces Section 3.4 (Mutation Spectra) of:
#   Dodge et al. (2025) "MutSeqR: an open source R package for standardized
#   analysis of error-corrected next-generation sequencing data in genetic
#   toxicology". Bioinformatics Advances, vbaf265.
#
# A mutation spectrum is the pattern of mutation subtypes within a group.
# Comparing spectra between BaP dose groups and controls reveals the
# characteristic mutational fingerprint of BaP exposure.
#
# TWO ANALYSES
# ------------
# 1. spectra_comparison(): tests whether the proportion of each mutation
#    subtype differs significantly between each BaP dose group and the
#    control. Uses a G2 likelihood ratio statistic (chi-squared for large
#    N, F-distribution when N/(R-1) < 20).
#
# 2. cluster_spectra(): hierarchically clusters samples by their subtype
#    proportions using cosine distance and Ward linkage. Samples from
#    the same dose group should cluster together if BaP induces a
#    dose-dependent shift in the mutation spectrum.
#
# DATA REQUIREMENTS
# ------------------
# Both functions require MF data at a subtype resolution. We use:
#   - "base_6": 6 SNV subtypes (C>A, C>G, C>T, T>A, T>C, T>G) in
#     pyrimidine reference, plus non-SNV variation types.
#   - Grouped by dose_group for spectra_comparison()
#   - Grouped by new_sample_id for cluster_spectra()
#
# NOTE ON DATA SOURCE
# --------------------
# The spectra_comparison() example in the docs uses eh[["EH9861"]]
# (filtered mutation data) rather than the processed EH9860 we used
# in Scripts 01-04. We load directly from EH9860 (which already has
# filter_mut applied) to maintain pipeline consistency.
#
# REPRODUCTION TARGETS (section 3.4, Figures 6A and 6B)
# -------------------------------------------------------
# - BaP subtype proportions significantly different from control at
#   ALL three dose groups (P < .001)
# - Dose-dependent increase in C:G > A:T mutation proportions
# - Samples cluster by dose group in hierarchical clustering
# - High and medium dose samples separate from control and low dose
#
# =============================================================================

library(MutSeqR)
library(ExperimentHub)
library(dplyr)

cat("=== Script 05: Mutation Spectra Analysis ===\n\n")

# ---- 1. Load filtered data --------------------------------------------------
# We use the filtered dataset saved in Script 01 to maintain consistency
# with our filtering decisions (2660 flagged rows).

bap_filtered <- readRDS("outputs/01_bap_filtered.rds")
cat("Loaded filtered data:", nrow(bap_filtered), "rows\n")
cat("Flagged rows:", sum(bap_filtered$filter_mut), "\n\n")

# ---- 2. Calculate MF at 6-base subtype resolution by dose group -------------
# spectra_comparison() requires MF data grouped by the experimental variable
# (dose_group) with subtype resolution. We use "base_6" (6 SNV subtypes in
# pyrimidine reference) — the resolution used in Figure 6A of the paper.
#
# MFMin is used (mf_type = "min") as the paper recommends for spectra
# analysis, because the independence assumption underlying the G2 test
# requires that mutations are counted as independent events.

cat("=== Calculating MF at base_6 resolution (by dose group) ===\n")
mf_dose_6 <- calculate_mf(
  mutation_data      = bap_filtered,
  cols_to_group      = "dose_group",
  subtype_resolution = "base_6",
  summary            = TRUE
)

cat("MF by dose group (base_6) dimensions:",
    nrow(mf_dose_6), "rows x", ncol(mf_dose_6), "columns\n")
cat("Subtypes present:\n")
print(table(mf_dose_6$normalized_subtype))
cat("\n")

# ---- 3. Calculate MF at 6-base resolution by sample (for clustering) --------
# cluster_spectra() operates on per-sample data so that individual animals
# can be clustered. We retain dose_group as metadata for interpretation.

cat("=== Calculating MF at base_6 resolution (by sample) ===\n")
mf_sample_6 <- calculate_mf(
  mutation_data        = bap_filtered,
  cols_to_group        = c("new_sample_id", "dose_group"),
  subtype_resolution   = "base_6",
  summary              = TRUE,
  retain_metadata_cols = "dose"
)

# Filter to real observations (same fix as Script 02)
mf_sample_6_real <- mf_sample_6 %>%
  filter(!is.na(dose))

cat("MF by sample (base_6) dimensions:",
    nrow(mf_sample_6_real), "rows (after filtering zeros)\n\n")

# ---- 4. Spectra comparison: each dose group vs control ---------------------
# We compare Low, Medium, and High vs Control.
# Col1 = treated group, Col2 = reference (control).

contrasts_spectra <- data.frame(
  col1 = c("Low", "Medium", "High"),
  col2 = rep("Control", 3)
)

cat("=== Running spectra_comparison() ===\n")
cat("Comparing each dose group vs Control at base_6 resolution\n\n")

spectra_results <- spectra_comparison(
  mf_data      = mf_dose_6,
  exp_variable = "dose_group",
  mf_type      = "min",
  contrasts     = contrasts_spectra
)

cat("=== SPECTRA COMPARISON RESULTS ===\n")
print(spectra_results)

# ---- 5. Reproduction check: significance at all doses ----------------------

cat("\n=== REPRODUCTION CHECK ===\n")
cat("Target: significant spectrum difference at ALL dose groups (P < .001)\n\n")

if (!is.null(spectra_results)) {
  for (i in 1:nrow(spectra_results)) {
    row <- spectra_results[i, ]
    sig <- ifelse(row$p.value < 0.001, "***",
                  ifelse(row$p.value < 0.01, "**",
                         ifelse(row$p.value < 0.05, "*", "ns")))
    cat(sprintf("  %s: G2 = %.2f, p = %.2e  %s\n",
                rownames(spectra_results)[i],
                row$G2, row$p.value, sig))
  }
  all_sig <- all(spectra_results$p.value < 0.001)
  cat(sprintf("\nAll comparisons P < .001: %s | Target: TRUE\n", all_sig))
}

# ---- 6. C:G > A:T dose-dependent increase ----------------------------------
# The paper reports dose-dependent increases in C:G > A:T proportions,
# consistent with BaP's known mechanism (bulky adducts at guanine N2).

cat("\n=== C:G > A:T PROPORTION BY DOSE GROUP ===\n")
cat("(normalized proportion of C>A mutations out of all SNVs)\n\n")

ca_proportions <- mf_dose_6 %>%
  filter(normalized_subtype == "C>A") %>%
  select(dose_group, proportion_min) %>%
  mutate(dose_group = factor(dose_group,
                             levels = c("Control", "Low", "Medium", "High"))) %>%
  arrange(dose_group)

print(ca_proportions)
cat("\nExpected trend: increasing C>A proportion with BaP dose\n")

# ---- 7. Hierarchical clustering of samples ---------------------------------
# cluster_spectra() groups samples by cosine distance of their 6-base
# subtype proportions. The paper (Figure 6B) shows samples clustering
# by dose group, with high/medium separating from control/low.

cat("\n=== Hierarchical Clustering ===\n")

# Filter to SNV subtypes for clustering (as per paper Figure 6B)
mf_snv_6 <- mf_sample_6_real %>%
  filter(normalized_subtype %in% c("C>A", "C>G", "C>T", "T>A", "T>C", "T>G"))

dendro <- cluster_spectra(
  mf_data      = mf_snv_6,
  group_col    = "new_sample_id",
  response_col = "proportion_min",
  subtype_col  = "normalized_subtype",
  dist         = "euclidean",
  cluster_method = "ward.D"
)

cat("Dendrogram created successfully\n")
cat("Cluster order (left to right):\n")
print(labels(dendro))

cat("\nExpected: High and Medium samples cluster separately from\n")
cat("Control and Low samples (consistent with Figure 6B)\n")

# ---- 8. Save results --------------------------------------------------------

dir.create("outputs", showWarnings = FALSE)
saveRDS(mf_dose_6,        "outputs/05_mf_dose_6.rds")
saveRDS(mf_sample_6_real, "outputs/05_mf_sample_6.rds")
saveRDS(spectra_results,  "outputs/05_spectra_results.rds")
saveRDS(dendro,           "outputs/05_dendrogram.rds")

if (!is.null(spectra_results)) {
  write.csv(spectra_results, "outputs/05_spectra_results.csv")
}
write.csv(ca_proportions,   "outputs/05_ca_proportions.csv",
          row.names = FALSE)

cat("\n=== Results saved ===\n")
cat("outputs/05_mf_dose_6.rds\n")
cat("outputs/05_mf_sample_6.rds\n")
cat("outputs/05_spectra_results.rds / .csv\n")
cat("outputs/05_dendrogram.rds\n")
cat("outputs/05_ca_proportions.csv\n")
cat("\nReady for Script 06: COSMIC Signature Fitting\n")

# =============================================================================
# END OF SCRIPT 05
# =============================================================================