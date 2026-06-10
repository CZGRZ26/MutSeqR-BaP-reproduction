# =============================================================================
# MutSeqR BaP Reproduction — Project 1
# Script 01: Data Import & Variant Filtering
# =============================================================================
#
# OVERVIEW
# --------
# This script reproduces Section 3.1 (Variant Filtering) of:
#   Dodge et al. (2025) "MutSeqR: an open source R package for standardized
#   analysis of error-corrected next-generation sequencing data in genetic
#   toxicology". Bioinformatics Advances, vbaf265.
#   https://doi.org/10.1093/bioadv/vbaf265
#
# The dataset is Project 1: MutaMouse male bone marrow samples exposed to
# benzo[a]pyrene (BaP) at four doses (Control, Low, Medium, High), sequenced
# using TwinStrand Duplex Sequencing on the 20-target Mouse Mutagenesis Panel
# (48 kb total). There are 24 samples across the 4 dose groups (6 per group).
#
# REPRODUCTION TARGETS (§3.1)
# ---------------------------
# The paper reports that filter_mut() flagged or removed 2,660 out of
# 1,152,911 rows of mutation data, broken down as:
#
#   ~612  putative germline variants (VAF > 0.01)
#     20  SNVs overlapping germline MNVs
#   2021  custom filter (EndRepairFillInArtifact)
#     22  records outside target regions
#   ----
#   2660  total unique flagged rows
#
# Note: individual counts sum to 2,675. The difference of 15 reflects
# rows satisfying more than one criterion simultaneously (double-counted
# when summed naively):
#   - 4 rows: germline AND outside target regions
#   - 11 rows: germline AND SNV-in-germline-MNV
#
# RESULT: All four individual counts and the total of 2,660 reproduced exactly.
#
# NOTE ON FILTER NAME DISCREPANCY
# --------------------------------
# The paper refers to the custom filter as "EndRepairFillinArtifact"
# (lowercase 'i' in 'In'). The actual data uses "EndRepairFillInArtifact"
# (capital 'I'). This is a minor discrepancy between the paper text and the
# data; grepl() is used rather than exact matching to handle this and to
# capture compound filter strings like "EndRepairFillInArtifact,NM8.0".
#
# =============================================================================

# ---- 0. Load required libraries ---------------------------------------------

library(MutSeqR)       # Main analysis package (v1.1.0, working_version branch)
library(ExperimentHub) # Bioconductor data retrieval
library(dplyr)         # Data manipulation

# ---- 1. Capture session info for reproducibility ----------------------------
# Recording the session info at the start ensures anyone re-running this
# script can identify any package version differences that might cause
# numeric discrepancies.

cat("=== Session Info ===\n")
cat("Date:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("R version:", R.version$version.string, "\n")
cat("MutSeqR version:", as.character(packageVersion("MutSeqR")), "\n\n")

# ---- 2. Load the BaP Project 1 data from ExperimentHub ---------------------
# MutSeqR ships example data via Bioconductor's ExperimentHub. The full
# processed Project 1 dataset (all 24 samples, all 1,152,911 rows) is stored
# under record EH9860 "Example mutation data".
#
# Note: EH9857 ("Example import_mut_data") is a single-sample subset only —
# it does NOT contain all 24 samples and cannot be used to reproduce the
# paper's filtering counts. EH9860 is the correct record for this analysis.
#
# EH9860 is a fully annotated tibble (38 columns) produced by running
# import_mut_data() on the raw tabular files from TwinStrand's pipeline.
# It includes metadata columns added by MutSeqR (dose_group, variation_type,
# vaf, in_regions, filter, etc.) that are needed for filtering.

cat("=== Loading Data ===\n")
eh <- ExperimentHub()
bap_raw <- eh[["EH9860"]]

cat("Dimensions:", dim(bap_raw)[1], "rows x", dim(bap_raw)[2], "columns\n")
cat("Unique samples:", length(unique(bap_raw$sample)), "\n")
cat("Dose groups:\n")
print(table(bap_raw$dose_group))
cat("\n")

# ---- 3. Apply variant filters -----------------------------------------------
# The paper applies four filters from §2.3.1, following TwinStrand's
# recommendations for Duplex Sequencing data. We reproduce each filter as a
# logical flag vector, then combine them with OR to get unique flagged rows.
# Using flags rather than removing rows preserves total_depth values for
# downstream mutation frequency calculations — this mirrors how MutSeqR's
# filter_mut() works internally.

cat("=== Applying Filters ===\n")

# Filter 1: Putative germline variants (VAF > 0.01)
# Variants with a variant allele fraction above 1% are likely inherited
# germline polymorphisms rather than somatic mutations induced by BaP.
# These would inflate mutation frequency estimates if included.
flag_germline <- bap_raw$vaf > 0.01
cat("1. Putative germline (VAF > 0.01):",
    sum(flag_germline, na.rm = TRUE), "| Target: ~612\n")

# Filter 2: EndRepairFillInArtifact (custom TwinStrand filter)
# During Duplex Sequencing library preparation, end-repair and fill-in steps
# can introduce false mutations at the ends of DNA fragments. TwinStrand's
# pipeline flags these in the 'filter' column. We use grepl() rather than
# exact matching because the filter column contains compound strings
# (e.g. "EndRepairFillInArtifact,NM8.0"). Note: the paper text uses the
# spelling "EndRepairFillinArtifact" (lowercase 'i') but the data uses
# "EndRepairFillInArtifact" (capital 'I').
flag_end_repair <- grepl("EndRepairFillInArtifact", bap_raw$filter)
cat("2. EndRepairFillInArtifact:       ",
    sum(flag_end_repair, na.rm = TRUE), "| Target: 2021\n")

# Filter 3: Outside target regions
# The Mouse Mutagenesis Panel targets 20 specific genomic regions (2.4 kb
# each). Records outside these regions are not part of the intended sequencing
# panel and must be excluded. MutSeqR annotates each row with an 'in_regions'
# logical column during import.
flag_outside <- !bap_raw$in_regions
cat("3. Outside target regions:        ",
    sum(flag_outside, na.rm = TRUE), "| Target: 22\n")

# Filter 4: SNVs overlapping germline MNVs
# Germline multi-nucleotide variants (MNVs) can create false sub-clonal SNVs
# during variant calling: no-calls in reads supporting the germline MNV may
# appear as minor haplotypes (apparent SNVs). We identify germline MNVs as
# MNVs with VAF > 0.01, then flag any SNV whose position falls within the
# span of one of those germline MNVs.
germline_mnvs <- bap_raw %>%
  filter(variation_type == "mnv", vaf > 0.01)

germ_mnv_pos <- germline_mnvs %>%
  select(contig, start, end) %>%
  distinct()

flag_snv_mnv <- bap_raw$variation_type == "snv" &
  paste(bap_raw$contig, bap_raw$start) %in%
  unlist(lapply(1:nrow(germ_mnv_pos), function(i) {
    paste(germ_mnv_pos$contig[i],
          seq(germ_mnv_pos$start[i], germ_mnv_pos$end[i]))
  }))
cat("4. SNVs in germline MNVs:         ",
    sum(flag_snv_mnv, na.rm = TRUE), "| Target: 20\n")

# ---- 4. Combine flags and count unique flagged rows -------------------------
# Summing the four individual counts gives 2,675 — 15 more than the paper's
# 2,660. The excess reflects rows satisfying more than one criterion:
#   - 4 rows: germline AND outside target regions
#   - 11 rows: germline AND SNV-in-germline-MNV
# Using OR (|) correctly counts each flagged row only once.

cat("\n=== Overlap Between Filters ===\n")
cat("germline & end_repair:  ", sum(flag_germline & flag_end_repair, na.rm = TRUE), "\n")
cat("germline & outside:     ", sum(flag_germline & flag_outside,    na.rm = TRUE), "\n")
cat("germline & snv_mnv:     ", sum(flag_germline & flag_snv_mnv,    na.rm = TRUE), "\n")
cat("end_repair & outside:   ", sum(flag_end_repair & flag_outside,  na.rm = TRUE), "\n")
cat("end_repair & snv_mnv:   ", sum(flag_end_repair & flag_snv_mnv,  na.rm = TRUE), "\n")
cat("outside & snv_mnv:      ", sum(flag_outside & flag_snv_mnv,     na.rm = TRUE), "\n")

total_unique <- sum(
  flag_germline | flag_end_repair | flag_outside | flag_snv_mnv,
  na.rm = TRUE
)

cat("\n=== Final Filtering Tally ===\n")
cat("Putative germline (VAF > 0.01): ",
    sum(flag_germline, na.rm = TRUE),   "| Target: ~612\n")
cat("EndRepairFillInArtifact:        ",
    sum(flag_end_repair, na.rm = TRUE), "| Target: 2021\n")
cat("Outside target regions:         ",
    sum(flag_outside, na.rm = TRUE),    "| Target: 22\n")
cat("SNVs in germline MNVs:          ",
    sum(flag_snv_mnv, na.rm = TRUE),    "| Target: 20\n")
cat("TOTAL (unique):                 ",
    total_unique,                       "| Target: 2660\n")

# ---- 5. Add filter flag to dataset and save ---------------------------------
# We add a single combined filter_mut column (TRUE = flagged, exclude from
# mutation counts) and save the result for use in Script 02 (mutation
# frequency calculation). Flagged rows are retained — not removed — so their
# depth values remain available for denominator calculations.

bap_filtered <- bap_raw %>%
  mutate(filter_mut = flag_germline | flag_end_repair |
           flag_outside   | flag_snv_mnv)

cat("\nRows flagged (filter_mut == TRUE):",
    sum(bap_filtered$filter_mut, na.rm = TRUE), "\n")
cat("Rows retained for analysis:      ",
    sum(!bap_filtered$filter_mut, na.rm = TRUE), "\n")

dir.create("outputs", showWarnings = FALSE)
saveRDS(bap_filtered, "outputs/01_bap_filtered.rds")
cat("\nFiltered dataset saved to: outputs/01_bap_filtered.rds\n")
cat("Ready for Script 02: Mutation Frequency Calculation\n")

# =============================================================================
# END OF SCRIPT 01
# =============================================================================