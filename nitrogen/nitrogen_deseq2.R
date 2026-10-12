#!/usr/bin/env Rscript
# Nitrogen functional-module expression analysis.
# Run: Rscript "/path/to/nitrogen/nitrogen_deseq2.R"
# Software: DESeq2 1.52.0 and base R; no packages are installed by this script.
# Required inputs in the same directory as this script:
#   all_gene_raw_count_matrix_for_DESeq2.tsv.gz
#   nitrogen_gene_total_count_wide.csv
#   nitrogen_gene_tpm.csv
# Both comparisons include paired biological replicates B1, B2, and B3:
# microaerobic S1 -> E1; anoxic E1 -> E2; n = 3 pairs per comparison.
# Outputs generated at runtime in results/:
#   Table_S7_numeric_all_rows.csv and sessionInfo.txt

options(stringsAsFactors = FALSE, width = 140)
if (!requireNamespace("DESeq2", quietly = TRUE)) {
  stop("DESeq2 1.52.0 is required. Install it before running this script.")
}
suppressPackageStartupMessages(library(DESeq2))
if (as.character(packageVersion("DESeq2")) != "1.52.0") {
  warning("This analysis specifies DESeq2 1.52.0; the installed version is ", packageVersion("DESeq2"), ".")
}
script_argument <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(script_argument) != 1L) stop("Run this file with Rscript.")
script_value <- sub("^--file=", "", script_argument)
# R can encode spaces as ~+~ in the --file argument.
if (!file.exists(script_value)) script_value <- gsub("~+~", " ", script_value, fixed = TRUE)
script_path <- normalizePath(script_value, mustWork = TRUE)
input_dir <- dirname(script_path)
output_dir <- file.path(input_dir, "results")
input_names <- c("all_gene_raw_count_matrix_for_DESeq2.tsv.gz",
  "nitrogen_gene_total_count_wide.csv", "nitrogen_gene_tpm.csv")
input_files <- file.path(input_dir, input_names)
if (!all(file.exists(input_files))) {
  stop("Missing input files: ", paste(input_names[!file.exists(input_files)], collapse = ", "))
}
connection <- gzfile(input_files[1], "rt")
all_gene_counts <- read.delim(connection, check.names = FALSE)
close(connection)
nitrogen_counts <- read.csv(input_files[2], check.names = FALSE)
gene_tpm <- read.csv(input_files[3], check.names = FALSE)

# Each row defines a gene-to-module assignment within a functional group.
# MAG denotes a functional group here; several groups comprise multiple MAGs.
module_map <- read.csv(text = "MAG,Gene,Module
AOB,amoA,amoABC
AOB,amoB,amoABC
AOB,amoC,amoABC
AOB,hao,hao
CMX,amoA,amoABC
CMX,amoB,amoABC
CMX,amoC,amoABC
CMX,hao,hao
CMX,nxrA,nxrAB
CMX,nxrB,nxrAB
NOB,nxrA,nxrAB
NOB,nxrB,nxrAB
AMX,hao,hao
AMX,hzsA,hzsABC
AMX,hzsB,hzsABC
AMX,hzsC,hzsABC
AMX,hzo,hzo
AMX,narG,narGH
AMX,narH,narGH
AMX,nrfA,nrfAH
AMX,nrfH,nrfAH
Complete_DNRA_MAGs,narG,narGH
Complete_DNRA_MAGs,narH,narGH
Complete_DNRA_MAGs,napA,napAB
Complete_DNRA_MAGs,napB,napAB
Complete_DNRA_MAGs,nrfA,nrfAH
Complete_DNRA_MAGs,nrfH,nrfAH
Complete_DNRA_MAGs,nirB,nirBD
Complete_DNRA_MAGs,nirD,nirBD
Complete_Denitrifiers_MAGs,narG,narGH
Complete_Denitrifiers_MAGs,narH,narGH
Complete_Denitrifiers_MAGs,napA,napAB
Complete_Denitrifiers_MAGs,napB,napAB
Complete_Denitrifiers_MAGs,nirK,nirKS
Complete_Denitrifiers_MAGs,nirS,nirKS
Complete_Denitrifiers_MAGs,norB,norBC
Complete_Denitrifiers_MAGs,norC,norBC
Complete_Denitrifiers_MAGs,nosZ,nosZ
Partial_denitrification_MAGs,narG,narGH
Partial_denitrification_MAGs,narH,narGH
Partial_denitrification_MAGs,napA,napAB
Partial_denitrification_MAGs,napB,napAB", check.names = FALSE)
stopifnot(!anyDuplicated(all_gene_counts$FeatureID),
  !anyDuplicated(module_map[c("MAG", "Gene")]))
normalize_group <- function(x) {
  x[x == "Complete_Denitrification_and_DNRA_MAGs"] <- "Complete_Denitrifiers_MAGs"
  x
}
normalize_gene <- function(x) {
  x[x == "nxrA_or_narG"] <- "nxrA"
  x[x == "nxrB_or_narH"] <- "nxrB"
  x
}
nitrogen_counts <- nitrogen_counts[
  !is.na(nitrogen_counts$target_gene) & nitrogen_counts$target_gene != "" &
    !is.na(nitrogen_counts$functional_group) & nitrogen_counts$functional_group != "", ]
nitrogen_counts$MAG <- normalize_group(nitrogen_counts$functional_group)
nitrogen_counts$Gene <- normalize_gene(nitrogen_counts$target_gene)
nitrogen_counts$FeatureID <- with(nitrogen_counts,
  paste(contig, start, end, strand, locus_tag, sep = "|"))
stopifnot(!anyDuplicated(nitrogen_counts$FeatureID))
gene_tpm <- gene_tpm[
  !is.na(gene_tpm$Gene) & gene_tpm$Gene != "" & !is.na(gene_tpm$MAG) & gene_tpm$MAG != "", ]
gene_tpm$MAG <- normalize_group(gene_tpm$MAG)
gene_tpm$Gene <- normalize_gene(gene_tpm$Gene)
samples_all <- grep("^AJP15_.*_RNA$", names(all_gene_counts), value = TRUE)
required_samples <- as.vector(outer(paste0("AJP15_B", 1:3),
  c("S1_RNA", "E1_RNA", "E2_RNA"), paste0))
for (object in list(all_gene_counts, nitrogen_counts, gene_tpm)) {
  stopifnot(all(required_samples %in% names(object)))
}
assert_integer_counts <- function(x, description) {
  m <- as.matrix(x)
  if (!is.numeric(m) || anyNA(m) || any(!is.finite(m)) || any(m < 0) || any(m != round(m))) {
    stop(description, " must contain finite, nonnegative integer counts without missing values.")
  }
}
assert_integer_counts(all_gene_counts[, samples_all, drop = FALSE], "All-gene matrix")
assert_integer_counts(nitrogen_counts[, samples_all, drop = FALSE], "Nitrogen-gene matrix")
reference_index <- match(nitrogen_counts$FeatureID, all_gene_counts$FeatureID)
stopifnot(!anyNA(reference_index))
stopifnot(all(as.matrix(nitrogen_counts[, samples_all, drop = FALSE]) ==
  as.matrix(all_gene_counts[reference_index, samples_all, drop = FALSE])))

# Sum raw counts over mapped loci, subunits and MAGs within each
# functional-group/module/sample combination. Biological samples are not pooled.
mapped_loci <- merge(nitrogen_counts, module_map, by = c("MAG", "Gene"), all = FALSE, sort = FALSE)
stopifnot(!anyDuplicated(mapped_loci$FeatureID))
module_ids <- paste(mapped_loci$MAG, mapped_loci$Module, sep = "__")
module_counts <- rowsum(round(as.matrix(mapped_loci[, samples_all, drop = FALSE])),
  group = module_ids, reorder = TRUE)
module_counts <- module_counts[order(rownames(module_counts), method = "radix"), , drop = FALSE]
storage.mode(module_counts) <- "integer"
module_annotation <- data.frame(FeatureID = rownames(module_counts),
  MAG = sub("__.*$", "", rownames(module_counts)),
  Module = sub("^.*__", "", rownames(module_counts)))
# Retain the computational identifier alongside the manuscript's displayed label.
module_annotation$Module_display <- ifelse(
  module_annotation$MAG == "AMX" & module_annotation$Module == "narGH", "nxrAB", module_annotation$Module)

# TPM is a separate descriptive quantity and is never input to DESeq2.
# The input TPM table contains functional-group/gene sums; module subunits
# are summed within each sample before computing the replicate mean, SD and SE.
tpm_mapped <- merge(gene_tpm, module_map, by = c("MAG", "Gene"), all = FALSE, sort = FALSE)
module_tpm <- rowsum(as.matrix(tpm_mapped[, samples_all, drop = FALSE]),
  group = paste(tpm_mapped$MAG, tpm_mapped$Module, sep = "__"), reorder = TRUE)
stopifnot(setequal(rownames(module_counts), rownames(module_tpm)))
module_tpm <- module_tpm[module_annotation$FeatureID, , drop = FALSE]
stopifnot(!anyNA(module_tpm), all(is.finite(module_tpm)), all(module_tpm >= 0))
significance_label <- function(p) {
  ifelse(is.na(p), NA_character_, ifelse(p < 0.001, "***",
    ifelse(p < 0.01, "**", ifelse(p < 0.05, "*", ""))))
}
run_stage <- function(stage_name, start_suffix, end_suffix) {
  biological_replicates <- paste0("B", 1:3)
  selected_samples <- as.vector(rbind(
    paste0("AJP15_", biological_replicates, start_suffix, "_RNA"),
    paste0("AJP15_", biological_replicates, end_suffix, "_RNA")))
  coldata <- data.frame(
    BioRep = factor(rep(biological_replicates, each = 2), levels = biological_replicates),
    Stage = factor(rep(c("Start", "End"), 3), levels = c("Start", "End")),
    row.names = selected_samples)
  count_all <- as.matrix(all_gene_counts[, selected_samples, drop = FALSE])
  rownames(count_all) <- all_gene_counts$FeatureID
  storage.mode(count_all) <- "integer"
  stopifnot(identical(rownames(coldata), colnames(count_all)))
  # Estimate size factors from all genes for the six libraries of this comparison.
  # Transfer those factors to the functional-module count dataset.
  dds_all <- DESeqDataSetFromMatrix(countData = count_all, colData = coldata,
    design = ~ BioRep + Stage)
  dds_all <- estimateSizeFactors(dds_all, type = "ratio")
  sf <- sizeFactors(dds_all)
  dds <- DESeqDataSetFromMatrix(countData = module_counts[, selected_samples, drop = FALSE],
    colData = coldata, design = ~ BioRep + Stage)
  sizeFactors(dds) <- sf[colnames(dds)]
  keep <- rowSums(counts(dds) >= 10) >= 3
  if (!any(keep)) stop("No nitrogen modules pass filtering in ", stage_name)
  # Wald testing uses a parametric dispersion trend and DESeq2's default
  # fallback and Cook's-distance handling. No log-fold-change shrinkage is used.
  dds <- DESeq(dds[keep, ])
  res <- results(dds, contrast = c("Stage", "End", "Start"),
    independentFiltering = FALSE, pAdjustMethod = "BH")
  result <- data.frame(FeatureID = rownames(res), as.data.frame(res), row.names = NULL)
  idx <- match(module_annotation$FeatureID, result$FeatureID)
  full_result <- cbind(module_annotation, Retained = keep,
    result[idx, setdiff(names(result), "FeatureID"), drop = FALSE])
  rownames(full_result) <- NULL
  full_result$Test_status <- ifelse(!keep, "Not_tested_low_count",
    ifelse(is.na(full_result$pvalue), "P_value_unavailable",
      ifelse(is.na(full_result$padj), "Adjusted_P_value_unavailable", "Tested")))
  full_result$Significant_label <- significance_label(full_result$padj)
  endpoint_samples <- paste0("AJP15_", biological_replicates, end_suffix, "_RNA")
  endpoint_values <- module_tpm[, endpoint_samples, drop = FALSE]
  endpoint <- data.frame(module_annotation, Stage = stage_name, N = 3,
    B1_endpoint_TPM = endpoint_values[, 1], B2_endpoint_TPM = endpoint_values[, 2],
    B3_endpoint_TPM = endpoint_values[, 3], Endpoint_mean_TPM = rowMeans(endpoint_values),
    Endpoint_SD_TPM = apply(endpoint_values, 1, sd), row.names = NULL)
  # Sample SD uses denominator n-1. SE is SD/sqrt(n), not the model's lfcSE.
  endpoint$Endpoint_SE_TPM <- endpoint$Endpoint_SD_TPM / sqrt(endpoint$N)
  combined <- cbind(endpoint, full_result[, c("Retained", "baseMean", "log2FoldChange", "lfcSE",
    "stat", "pvalue", "padj", "Test_status", "Significant_label"), drop = FALSE])
  message(stage_name, ": ", sum(keep), "/", length(keep), " modules tested; n = 3 pairs; dispersion fit = ",
    attr(dispersionFunction(dds), "fitType"))
  combined
}

microaerobic <- run_stage("Microaerobic", "S1", "E1")
anoxic <- run_stage("Anoxic", "E1", "E2")
complete_results <- rbind(microaerobic, anoxic)
table_order <- c("AOB__amoABC", "AOB__hao", "CMX__amoABC", "CMX__hao", "CMX__nxrAB",
  "NOB__nxrAB", "AMX__hao", "AMX__hzsABC", "AMX__hzo", "AMX__narGH", "AMX__nrfAH",
  "Complete_DNRA_MAGs__narGH", "Complete_DNRA_MAGs__napAB", "Complete_DNRA_MAGs__nrfAH",
  "Complete_DNRA_MAGs__nirBD", "Complete_Denitrifiers_MAGs__narGH", "Complete_Denitrifiers_MAGs__napAB",
  "Complete_Denitrifiers_MAGs__nirKS", "Complete_Denitrifiers_MAGs__norBC", "Complete_Denitrifiers_MAGs__nosZ",
  "Partial_denitrification_MAGs__narGH", "Partial_denitrification_MAGs__napAB")
stopifnot(setequal(table_order, module_annotation$FeatureID))
complete_results <- complete_results[order(match(complete_results$Stage, c("Microaerobic", "Anoxic")),
  match(complete_results$FeatureID, table_order)), ]
rownames(complete_results) <- NULL
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
write.csv(complete_results, file.path(output_dir, "Table_S7_numeric_all_rows.csv"),
  row.names = FALSE, na = "NA")
capture.output(sessionInfo(), file = file.path(output_dir, "sessionInfo.txt"))
message("Completed nitrogen analysis.")
