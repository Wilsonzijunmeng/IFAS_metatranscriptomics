#!/usr/bin/env Rscript
# Carbon module analysis: three paired biological replicates (B1, B2, B3).
# Run: Rscript carbon_deseq2.R [input_directory]
# With no argument, inputs and results are located beside this script.
# Validated with R 4.6.1 and DESeq2 1.52.0; no packages are installed by this script.
# Three inputs:
#   acetate_propionate_gene_TPM_wide_main_only.tsv: locus annotations and TPM.
#   carbon_gene_raw_counts.csv: target-locus raw counts for module aggregation.
#   all_gene_raw_count_matrix_for_DESeq2.tsv.gz: all-gene normalization counts.
# Outputs: results/Table_S9_full_precision.csv and results/sessionInfo.txt.

args <- gsub("~+~", " ", commandArgs(trailingOnly = TRUE), fixed = TRUE)
if (length(args) > 1L) stop("Usage: Rscript carbon_deseq2.R [input_directory]")
if (length(args) == 1L) {
  input_dir <- normalizePath(args[1L], mustWork = TRUE)
} else {
  script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  script_arg <- gsub("~+~", " ", script_arg, fixed = TRUE)
  if (length(script_arg) != 1L) stop("Use Rscript or supply the input directory.")
  input_dir <- dirname(normalizePath(sub("^--file=", "", script_arg), mustWork = TRUE))
}
if (!requireNamespace("DESeq2", quietly = TRUE)) {
  stop("DESeq2 is required. This analysis was validated with DESeq2 1.52.0.")
}
suppressPackageStartupMessages(library(DESeq2))

input_names <- c("acetate_propionate_gene_TPM_wide_main_only.tsv",
                 "carbon_gene_raw_counts.csv",
                 "all_gene_raw_count_matrix_for_DESeq2.tsv.gz")
input_files <- file.path(input_dir, input_names)
if (!all(file.exists(input_files))) stop("All three carbon input files are required.")
input_hashes <- unname(tools::md5sum(input_files))
tpm <- read.delim(input_files[1L], check.names = FALSE, stringsAsFactors = FALSE)
target <- read.csv(input_files[2L], check.names = FALSE, stringsAsFactors = FALSE)
con <- gzfile(input_files[3L], open = "rt")
all_counts <- tryCatch(read.delim(con, check.names = FALSE, stringsAsFactors = FALSE),
                       finally = close(con))

# Map the nine AMX loci to seven modules; sum porA, porB and porC within samples.
module_order <- c("acs", "acdA", "acdB", "ackA", "pta", "xfp", "porABC")
module_map <- c(acs = "acs", acdA = "acdA", acdB = "acdB", ackA = "ackA",
                pta = "pta", xfp = "xfp", porA = "porABC", porB = "porABC", porC = "porABC")
required_tpm_columns <- c("Gene", "MAG", "contig", "start", "end", "strand", "locus_tag")
stopifnot(all(required_tpm_columns %in% names(tpm)),
          all(c("Gene", "Module", "FeatureID") %in% names(target)),
          "FeatureID" %in% names(all_counts))
tpm <- tpm[!is.na(tpm$Gene) & tpm$Gene != "" & !is.na(tpm$MAG) & tpm$MAG == "AMX", ]
tpm$Gene <- trimws(as.character(tpm$Gene))
tpm$FeatureID <- with(tpm, paste(contig, start, end, strand, locus_tag, sep = "|"))
tpm$Module <- unname(module_map[tpm$Gene])
stopifnot(nrow(tpm) == 9L, nrow(target) == 9L,
          !anyNA(tpm$Module), setequal(tpm$Gene, names(module_map)),
          !anyNA(tpm$FeatureID), !anyDuplicated(tpm$FeatureID),
          !anyNA(target$FeatureID), !anyDuplicated(target$FeatureID),
          !anyNA(all_counts$FeatureID), !anyDuplicated(all_counts$FeatureID))
target_index <- match(tpm$FeatureID, target$FeatureID)
all_index <- match(tpm$FeatureID, all_counts$FeatureID)
stopifnot(!anyNA(target_index), !anyNA(all_index))
target <- target[target_index, , drop = FALSE]
stopifnot(identical(target$Gene, tpm$Gene), identical(target$Module, tpm$Module))

# Check every target count against its exact FeatureID in the full matrix.
count_samples <- setdiff(names(all_counts), "FeatureID")
target_samples <- setdiff(names(target), c("Gene", "Module", "FeatureID"))
stopifnot(setequal(target_samples, count_samples))
target_counts <- as.matrix(target[, count_samples, drop = FALSE])
matched_all_counts <- as.matrix(all_counts[all_index, count_samples, drop = FALSE])
stopifnot(is.numeric(target_counts), !anyNA(target_counts),
          all(is.finite(target_counts)), all(target_counts >= 0),
          all(target_counts == round(target_counts)),
          identical(dim(target_counts), dim(matched_all_counts)),
          all(target_counts == matched_all_counts))

reps <- paste0("B", 1:3)
required_samples <- as.vector(sapply(reps, function(b) paste0("AJP15_", b, c("S1", "E1", "E2"), "_RNA")))
stopifnot(all(required_samples %in% count_samples), all(required_samples %in% names(tpm)))
tpm_values <- as.matrix(tpm[, required_samples, drop = FALSE])
stopifnot(is.numeric(tpm_values), all(is.finite(tpm_values)), all(tpm_values >= 0))
module_counts <- rowsum(target_counts, tpm$Module, reorder = FALSE)[module_order, , drop = FALSE]
module_tpm <- rowsum(tpm_values, tpm$Module, reorder = FALSE)[module_order, , drop = FALSE]

run_contrast <- function(start_timepoint, end_timepoint, stage_name) {
  samples <- as.vector(rbind(paste0("AJP15_", reps, start_timepoint, "_RNA"),
                            paste0("AJP15_", reps, end_timepoint, "_RNA")))
  coldata <- data.frame(BioRep = factor(rep(reps, each = 2L)),
                        Stage = factor(rep(c("Start", "End"), 3L), levels = c("Start", "End")),
                        row.names = samples)

  # Estimate median-of-ratios size factors using all genes in these six samples.
  full_matrix <- as.matrix(all_counts[, samples, drop = FALSE])
  rownames(full_matrix) <- all_counts$FeatureID
  full_matrix[is.na(full_matrix)] <- 0
  stopifnot(is.numeric(full_matrix), all(is.finite(full_matrix)),
            all(full_matrix >= 0), all(full_matrix == round(full_matrix)))
  storage.mode(full_matrix) <- "integer"
  stopifnot(identical(colnames(full_matrix), rownames(coldata)))
  dds_all <- DESeq2::DESeqDataSetFromMatrix(full_matrix, coldata, ~ BioRep + Stage)
  dds_all <- DESeq2::estimateSizeFactors(dds_all, type = "ratio")
  sf <- DESeq2::sizeFactors(dds_all)

  counts <- round(module_counts[, samples, drop = FALSE])
  storage.mode(counts) <- "integer"
  dds <- DESeq2::DESeqDataSetFromMatrix(counts, coldata, ~ BioRep + Stage)
  DESeq2::sizeFactors(dds) <- sf[colnames(dds)]
  keep <- rowSums(DESeq2::counts(dds) >= 10L) >= 3L
  if (!any(keep)) stop("No carbon modules meet the count filter.")
  message(stage_name, ": ", sum(keep), "/", length(keep), " modules pass >=10 counts in >=3 of 6 samples.")

  # Paired negative-binomial Wald test; no fold-change shrinkage.
  # BH adjustment is per contrast; retain default Cook's outlier handling.
  dds <- DESeq2::DESeq(dds[keep, ], test = "Wald", fitType = "parametric", betaPrior = FALSE)
  res <- as.data.frame(DESeq2::results(dds, contrast = c("Stage", "End", "Start"),
                                      independentFiltering = FALSE, pAdjustMethod = "BH"))
  res <- res[match(module_order, rownames(res)), , drop = FALSE]

  # Calculate module TPM summaries from three endpoint biological replicates.
  endpoint <- module_tpm[, paste0("AJP15_", reps, end_timepoint, "_RNA"), drop = FALSE]
  endpoint_sd <- apply(endpoint, 1L, stats::sd)
  p <- res$padj
  label <- ifelse(is.na(p), "", ifelse(p < 0.001, "***", ifelse(p < 0.01, "**", ifelse(p < 0.05, "*", ""))))
  data.frame(Module = module_order, Stage = stage_name, n = ncol(endpoint),
             Endpoint_mean_TPM = rowMeans(endpoint), Endpoint_SD = endpoint_sd,
             Endpoint_SE = endpoint_sd / sqrt(ncol(endpoint)),
             DESeq2_log2FC = res$log2FoldChange,
             DESeq2_fold_change = 2^res$log2FoldChange,
             Adjusted_p_value = p, Significant_label = label, row.names = NULL)
}

# Microaerobic: S1 -> E1. Anoxic: E1 -> E2. B1/B2/B3 are included in both.
table_s9 <- rbind(run_contrast("S1", "E1", "Microaerobic"),
                  run_contrast("E1", "E2", "Anoxic"))
stopifnot(identical(input_hashes, unname(tools::md5sum(input_files))))
out_dir <- file.path(input_dir, "results")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
write.csv(table_s9, file.path(out_dir, "Table_S9_full_precision.csv"), row.names = FALSE, na = "NA")
capture.output(sessionInfo(), file = file.path(out_dir, "sessionInfo.txt"))
cat("Completed: results/Table_S9_full_precision.csv and results/sessionInfo.txt\n")
