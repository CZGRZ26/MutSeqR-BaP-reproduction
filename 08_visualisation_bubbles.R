# =============================================================================
# MutSeqR BaP Reproduction — Project 1
# Script 08: Bubble Plot Visualisation (§3.7)
# =============================================================================
#
# OVERVIEW
# --------
# This script reproduces Section 3.7 (Visualization of Multiplets MFMax) of:
#   Dodge et al. (2025) "MutSeqR: an open source R package for standardized
#   analysis of error-corrected next-generation sequencing data in genetic
#   toxicology". Bioinformatics Advances, vbaf265.
#
# BACKGROUND: MFMin vs MFMax
# ---------------------------
# MFMin counts each unique mutation position once per sample — this is the
# biologically conservative estimate assuming multiple reads with the same
# mutation reflect clonal expansion of a single event.
#
# MFMax counts every read supporting a non-reference allele — this is the
# upper bound and inflates when individual mutations are supported by many
# reads (multiplets). Large MFMax/MFMin ratios indicate clonal expansion.
#
# In this dataset, MFMax is considerably higher than MFMin across all dose
# groups, raising the question: is the inflation driven by a few very large
# multiplets, or by many moderately recurrent mutations?
#
# BUBBLE PLOT INTERPRETATION
# ---------------------------
# plot_bubbles() represents each mutation as a circle scaled by alt_depth
# (number of reads supporting the mutation). Large circles = highly recurrent
# multiplets. The distribution of circle sizes answers the question above.
#
# PAPER FINDINGS (§3.7, Figure 8)
# ---------------------------------
# Control: slight MFMin/MFMax difference driven by THREE large multiplets
#   - Two C:G>T:A mutations
#   - One deletion
#
# BaP dose groups: multiplets evenly distributed — no single mutation
# driving MFMax inflation. The increase in MFMax at BaP doses reflects
# many moderately recurrent mutations rather than a few large ones.
#
# =============================================================================
install.packages("packcircles")
library(packcircles)
library(MutSeqR)
library(dplyr)
library(ggplot2)

cat("=== Script 08: Bubble Plot Visualisation ===\n\n")

# ---- 1. Load filtered data --------------------------------------------------

bap_filtered <- readRDS("outputs/01_bap_filtered.rds")
cat("Loaded filtered data:", nrow(bap_filtered), "rows\n")
cat("filter_mut flagged:", sum(bap_filtered$filter_mut), "rows (auto-excluded by plot_bubbles)\n\n")

# ---- 2. Summarise MFMin vs MFMax by dose group ------------------------------
# Confirm the MFMax inflation observed in the paper

mf_by_sample <- readRDS("outputs/02_mf_by_sample.rds")

cat("=== MFMin vs MFMax ratio by dose group ===\n")
mf_by_sample %>%
  group_by(dose_group) %>%
  summarise(
    mean_mf_min = mean(mf_min_scaled),
    mean_mf_max = mean(mf_max_scaled),
    ratio       = mean(mf_max_scaled) / mean(mf_min_scaled),
    .groups = "drop"
  ) %>%
  mutate(dose_group = factor(dose_group,
                             levels = c("Control", "Low", "Medium", "High"))) %>%
  arrange(dose_group) %>%
  print()

# ---- 3. Identify the large multiplets in the control ------------------------
# Paper: three large multiplets in control — two C:G>T:A and one deletion

cat("\n=== Large multiplets in Control group ===\n")
control_multiplets <- bap_filtered %>%
  filter(dose_group == "Control",
         !filter_mut,
         alt_depth > 1) %>%
  arrange(desc(alt_depth)) %>%
  select(new_sample_id, contig, start, ref, alt, alt_depth,
         variation_type, normalized_subtype) %>%
  head(20)

cat("Top multiplets by alt_depth (Control):\n")
print(control_multiplets)

cat("\nLarge multiplets (alt_depth > 5) in Control:\n")
large_control <- bap_filtered %>%
  filter(dose_group == "Control", !filter_mut, alt_depth > 5)
print(large_control %>% select(new_sample_id, contig, start, ref, alt,
                               alt_depth, variation_type, normalized_subtype))

# Check paper's claim: 2 C:G>T:A + 1 deletion
cat("\nC>T mutations with high alt_depth in Control:\n")
ct_control <- bap_filtered %>%
  filter(dose_group == "Control", !filter_mut,
         normalized_subtype == "C>T", alt_depth > 1) %>%
  arrange(desc(alt_depth))
print(head(ct_control %>% select(new_sample_id, contig, start, ref, alt,
                                 alt_depth, normalized_subtype), 10))

cat("\nDeletions with high alt_depth in Control:\n")
del_control <- bap_filtered %>%
  filter(dose_group == "Control", !filter_mut,
         variation_type == "deletion", alt_depth > 1) %>%
  arrange(desc(alt_depth))
print(head(del_control %>% select(new_sample_id, contig, start, ref, alt,
                                  alt_depth, variation_type), 10))

# ---- 4. Generate bubble plots -----------------------------------------------
# Facet by dose_group to compare multiplet distributions across doses.
# color_by = "normalized_subtype" shows mutation type distribution.

cat("\n=== Generating bubble plots ===\n")

# Plot faceted by dose group
bubble_plot <- plot_bubbles(
  mutation_data  = bap_filtered,
  size_by        = "alt_depth",
  facet_col      = "dose_group",
  color_by       = "normalized_subtype",
  circle_spacing = 1
)

# Save plot
dir.create("outputs/figures", showWarnings = FALSE, recursive = TRUE)
ggsave(
  filename = "outputs/figures/08_bubble_plot_by_dose.png",
  plot     = bubble_plot,
  width    = 14,
  height   = 8,
  dpi      = 150
)

cat("Bubble plot saved to outputs/figures/08_bubble_plot_by_dose.png\n\n")

# ---- 5. Reproduction check --------------------------------------------------

cat("=== REPRODUCTION CHECK ===\n\n")

cat("Target: Control MFMax inflation driven by ~3 large multiplets\n")
cat("  (2 C:G>T:A mutations, 1 deletion)\n\n")

n_large_control <- nrow(bap_filtered %>%
                          filter(dose_group == "Control", !filter_mut, alt_depth > 5))
cat(sprintf("Large multiplets (alt_depth > 5) in Control: %d | Target: ~3\n",
            n_large_control))

cat("\nTarget: BaP dose groups show evenly distributed multiplets\n")
cat("  (no single mutation driving MFMax)\n\n")

# Compare max alt_depth across dose groups
max_alt_by_dose <- bap_filtered %>%
  filter(!filter_mut) %>%
  group_by(dose_group) %>%
  summarise(
    max_alt_depth  = max(alt_depth),
    mean_alt_depth = mean(alt_depth[alt_depth > 1]),
    n_multiplets   = sum(alt_depth > 1),
    .groups = "drop"
  ) %>%
  mutate(dose_group = factor(dose_group,
                             levels = c("Control", "Low", "Medium", "High"))) %>%
  arrange(dose_group)

cat("Alt_depth summary by dose group:\n")
print(max_alt_by_dose)

# ---- 6. Save results --------------------------------------------------------

ggsave(
  filename = "outputs/figures/08_bubble_plot_by_dose.pdf",
  plot     = bubble_plot,
  width    = 14,
  height   = 8
)

write.csv(control_multiplets,
          "outputs/08_control_large_multiplets.csv",
          row.names = FALSE)
write.csv(max_alt_by_dose,
          "outputs/08_alt_depth_summary.csv",
          row.names = FALSE)

cat("\n=== Results saved ===\n")
cat("outputs/figures/08_bubble_plot_by_dose.png\n")
cat("outputs/figures/08_bubble_plot_by_dose.pdf\n")
cat("outputs/08_control_large_multiplets.csv\n")
cat("outputs/08_alt_depth_summary.csv\n")

# =============================================================================
# END OF SCRIPT 08
# =============================================================================