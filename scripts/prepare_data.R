#! /usr/bin/env Rscript

# Enable error reporting and stack traces
options(error = function() { traceback(2); q(status = 1) })
options(warn = 1)  # print warnings as they occur

library(argparse)
library(tidyverse)
library(mia)
library(TreeSummarizedExperiment)

# Create a parser
parser <- ArgumentParser(description = "Prepare TSE objects from QIIME2 output files")

# Add arguments
parser$add_argument("--features", type = "character", required = TRUE, help = "Path to feature table QZA file (e.g., final_table.qza)")
parser$add_argument("--taxonomy", type = "character", required = TRUE, help = "Path to taxonomy QZA file (e.g., final_taxonomy.qza)")
parser$add_argument("--tree", type = "character", required = TRUE, help = "Path to rooted tree QZA file (e.g., rooted_tree.qza)")
parser$add_argument("--metadata", type = "character", required = TRUE, help = "Path to metadata TSV file")
parser$add_argument("--exclude_samples", type = "character", required = FALSE, default = NULL, help = "Optional file of sample IDs to drop before analysis, one per line (# starts a comment)")
parser$add_argument("--seqs", type = "character", required = TRUE, help = "Path to sequences QZA file (e.g., final_seqs.qza)")
parser$add_argument("--output_dir", type = "character", required = FALSE, default = ".", help = "Output directory for RDS files (default: current directory)")

# Create output directory if it doesn't exist
if (!dir.exists(parser$parse_args()$output_dir)) {
  dir.create(parser$parse_args()$output_dir, recursive = TRUE)
}

# Parse arguments
args <- parser$parse_args()
cat("✓ Parsed command line arguments\n")

# Print argument values for debugging
cat("Arguments:\n")
cat("  features:", args$features, "\n")
cat("  taxonomy:", args$taxonomy, "\n")
cat("  tree:", args$tree, "\n")
cat("  metadata:", args$metadata, "\n")
cat("  seqs:", args$seqs, "\n")
cat("  exclude_samples:", args$exclude_samples, "\n")
cat("  output_dir:", args$output_dir, "\n")

# Check if input files exist
input_files <- list(
  features = args$features,
  taxonomy = args$taxonomy,
  tree = args$tree,
  metadata = args$metadata,
  seqs = args$seqs,
  exclude_samples = args$exclude_samples
)
input_files <- input_files[!vapply(input_files, is.null, logical(1))]

for (name in names(input_files)) {
  if (!file.exists(input_files[[name]])) {
    stop("Input file not found: ", name, " = ", input_files[[name]])
  } else {
    cat("✓ Found", name, "file:", input_files[[name]], "\n")
  }
}

cat("=== Loading QIIME2 data ===\n")
tse <- tryCatch({
  mia::importQIIME2(
    assay.file = args$features,
    row.file = args$taxonomy,
    tree.file = args$tree,
    col.file = args$metadata,
    refseq.file = args$seqs,
    prefix.rm = TRUE
  )
},
 error = function(e) {
  cat("ERROR in importQIIME2: ", e$message, "\n")
  stop("Failed to import QIIME2 data: ", e$message)
})

# Factorization
rownames(colData(tse)) <- colData(tse)$SampleID

## Set factors
tse$Site <- factor(tse$Gut_Site, levels = c("Ileum", "Proximal", "Distal", "Excrement", "Lesion"), labels = c("Ileum", "Proximal", "Distal", "Feces", "Lesion"))
tse$Location <- factor(tse$Site, levels = c("Ileum", "Proximal", "Distal", "Feces", "Lesion"), labels = c("Gut", "Gut", "Gut", "Feces", "Gut"))
tse$Gut_Site <- NULL

tse$Group <- factor(tse$Group, levels = c("K", "A1", "A2", "C"))

tse$Age <- as.numeric(tse$Age)
tse$DNA_concentration <- as.numeric(tse$DNA_concentration)

## Merge_Batch: the per-run table the sample was QIIME2-processed and decontaminated in
## (1 = batch1.tsv, 2 = batch2.tsv, 3 = batch3.tsv).
tse$Merge_Batch <- factor(tse$Merge_Batch)

tse$Diagnosis <- factor(tse$Group, levels = c("K", "A1", "A2", "C"), labels = c("Control", "Adenoma", "Advanced adenoma", "Carcinoma"))

### ---Filtering---

# Filter Eukaryota and NA, (NA) first, because needed for TRUE/FALSE
tse <- tse[!is.na(rowData(tse)$domain), ]
tse <- tse[rowData(tse)$domain != "Eukaryota", ]

# Samples sequenced and deposited but left out of the analysis (low quality / preparation
# issues); listed in a file so no sample ID is hard-coded here.
if (!is.null(args$exclude_samples)) {
  samples_to_remove <- readLines(args$exclude_samples)
  samples_to_remove <- trimws(sub("#.*", "", samples_to_remove))
  samples_to_remove <- samples_to_remove[samples_to_remove != ""]
  unknown <- setdiff(samples_to_remove, colnames(tse))
  if (length(unknown) > 0) {
    stop("exclude_samples lists IDs that are not in the data: ", paste(unknown, collapse = ", "))
  }
  tse <- tse[, !colnames(tse) %in% samples_to_remove]
  cat("✓ Excluded", length(samples_to_remove), "samples;", ncol(tse), "remain\n")
}

# Fixed sample and feature order. The merged table lists samples in the order the batches
# reached the merge step, which changes between runs. The seeded rarefaction in r_prepare
# draws sample by sample, so without a fixed order the same seed gives different values.
tse <- tse[sort(rownames(tse)), sort(colnames(tse))]

# Patient is an identifier, not a number; order the levels naturally (2 before 10)
tse$Patient <- factor(tse$Patient, levels = stringr::str_sort(unique(tse$Patient), numeric = TRUE))

### --- Transform and prepare
tses <- list()

tses$full <- tse
tses$feces <- tse[, tse$Location %in% "Feces"]
tses$gut <- tse[, tse$Location %in% "Gut"]

## Reorder Site factor levels
tses$gut$Site <- factor(tses$gut$Site, levels = c("Ileum", "Proximal", "Distal", "Lesion"))
tses$full$Site <- factor(tses$full$Site, levels = c("Ileum", "Proximal", "Distal", "Lesion", "Feces"))

## Relative abundances (Bray-Curtis in r_prepare and the figures' taxon summaries use them)
tses <- lapply(tses, transformAssay, method = "relabundance")

saveRDS(tses$full, file = file.path(args$output_dir, "tse_full.rds"))
saveRDS(tses$feces, file = file.path(args$output_dir, "tse_feces.rds"))
saveRDS(tses$gut, file = file.path(args$output_dir, "tse_gut.rds"))
