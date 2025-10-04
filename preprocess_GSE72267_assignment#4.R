####################################################################
# preprocess_GSE72267_assignment.R
# Full preprocessing workflow for Affymetrix dataset (GSE72267)
# - QC before & after normalization (arrayQualityMetrics)
# - RMA normalization
# - Filtering low-intensity probes
# - Derive PD vs Control groups from phenotype
# - Differential expression (limma)
# - Saves CSVs, PNGs, and an assignment_summary.txt
####################################################################

## ---- Install/load ----
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
# uncomment to install missing packages on first run
BiocManager::install(c("GEOquery","affy","arrayQualityMetrics","limma","matrixStats"), ask=FALSE)
install.packages(c("dplyr","ggplot2"))

library(GEOquery)
library(affy)
library(arrayQualityMetrics)
library(matrixStats)
library(dplyr)
library(ggplot2)
library(limma)

## ---- Paths ----
workdir <- getwd()
raw_dir <- file.path(workdir, "Raw_GSE72267", "CEL_Files")
results_dir <- file.path(workdir, "Results")
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)

## ---- Read CEL files ----
cel_files <- list.celfiles(raw_dir, full.names = TRUE)
if (length(cel_files) == 0) stop("No CEL files found in: ", raw_dir, "\nExtract GSE72267_RAW.tar into this folder first.")
cat("Found", length(cel_files), "CEL files.\n")
affy_raw <- ReadAffy(filenames = cel_files)
print(affy_raw)

## ---- Programmatic QC (raw) ----
raw_probe_matrix_approx <- exprs(rma(affy_raw, normalize = FALSE, background = FALSE))
array_medians_raw <- apply(raw_probe_matrix_approx, 2, median, na.rm = TRUE)
q_raw <- quantile(array_medians_raw, c(0.25, 0.75)); iqr_raw <- IQR(array_medians_raw)
low_raw <- q_raw[1] - 1.5 * iqr_raw; high_raw <- q_raw[2] + 1.5 * iqr_raw
flag_med_raw <- (array_medians_raw < low_raw) | (array_medians_raw > high_raw)
pr_raw <- prcomp(t(raw_probe_matrix_approx), center = TRUE, scale. = TRUE)
scores_raw <- pr_raw$x[,1:5]
md_raw <- mahalanobis(scores_raw, colMeans(scores_raw), cov(scores_raw))
flag_pca_raw <- md_raw > (mean(md_raw) + 3 * sd(md_raw))
flag_raw_combined <- flag_med_raw | flag_pca_raw
raw_outlier_count_auto <- sum(flag_raw_combined)
raw_flag_df <- data.frame(sample = sampleNames(affy_raw), median = array_medians_raw,
                          mdist = md_raw, flag_median = flag_med_raw, flag_pca = flag_pca_raw,
                          flag_combined = flag_raw_combined, stringsAsFactors = FALSE)
write.csv(raw_flag_df, file = file.path(results_dir, "raw_array_flags_auto.csv"), row.names = FALSE)

## ---- arrayQualityMetrics raw (visual) ----
aqm_raw_dir <- file.path(results_dir, "arrayQM_raw"); dir.create(aqm_raw_dir, showWarnings = FALSE)
arrayQualityMetrics(expressionset = affy_raw, outdir = aqm_raw_dir, force = TRUE, do.logtransform = TRUE)
cat("arrayQualityMetrics (raw) saved to:", aqm_raw_dir, "\n")

## ---- RMA normalization ----
norm_affy <- rma(affy_raw)
norm_expr <- exprs(norm_affy)
cat("RMA complete. Normalized dims:", dim(norm_expr), "\n")

## ---- Programmatic QC (normalized) ----
norm_medians <- apply(norm_expr, 2, median, na.rm = TRUE)
q <- quantile(norm_medians, c(0.25, 0.75)); iqr <- IQR(norm_medians)
low_thr <- q[1] - 1.5 * iqr; high_thr <- q[2] + 1.5 * iqr
flag_med_norm <- (norm_medians < low_thr) | (norm_medians > high_thr)
pr_norm <- prcomp(t(norm_expr), center = TRUE, scale. = TRUE)
scores_norm <- pr_norm$x[,1:5]
mdist_norm <- mahalanobis(scores_norm, colMeans(scores_norm), cov(scores_norm))
flag_pca_norm <- mdist_norm > (mean(mdist_norm) + 3 * sd(mdist_norm))
flag_norm_combined <- flag_med_norm | flag_pca_norm
norm_outlier_count_auto <- sum(flag_norm_combined)
norm_flag_df <- data.frame(sample = colnames(norm_expr), median = norm_medians,
                           mdist = mdist_norm, flag_median = flag_med_norm, flag_pca = flag_pca_norm,
                           flag_combined = flag_norm_combined, stringsAsFactors = FALSE)
write.csv(norm_flag_df, file = file.path(results_dir, "norm_array_flags_auto.csv"), row.names = FALSE)

## ---- arrayQualityMetrics normalized (visual) ----
aqm_norm_dir <- file.path(results_dir, "arrayQM_normalized"); dir.create(aqm_norm_dir, showWarnings = FALSE)
arrayQualityMetrics(expressionset = norm_affy, outdir = aqm_norm_dir, force = TRUE, do.logtransform = FALSE)
cat("arrayQualityMetrics (normalized) saved to:", aqm_norm_dir, "\n")

## ---- Save boxplots & PCA ----
png(file.path(results_dir, "boxplot_raw.png"), width = 1400, height = 600)
boxplot(raw_probe_matrix_approx, main = "Raw intensities (approx)", las = 2)
dev.off()
png(file.path(results_dir, "boxplot_normalized.png"), width = 1400, height = 600)
boxplot(norm_expr, main = "RMA normalized", las = 2)
dev.off()

plot_and_save_pca <- function(mat, fname) {
  p <- prcomp(t(mat), center = TRUE, scale. = TRUE)
  df <- data.frame(PC1 = p$x[,1], PC2 = p$x[,2], Sample = rownames(p$x))
  png(file.path(results_dir, fname), width = 900, height = 700)
  plot(df$PC1, df$PC2, xlab = "PC1", ylab = "PC2", main = fname)
  text(df$PC1, df$PC2, labels = df$Sample, cex = 0.6, pos = 3)
  dev.off()
}
plot_and_save_pca(raw_probe_matrix_approx, "PCA_raw.png")
plot_and_save_pca(norm_expr, "PCA_normalized.png")

## ---- Filtering (expression-count method) ----
expr_thr <- log2(100)
keep_idx_expr <- rowSums(norm_expr > expr_thr) >= (0.25 * ncol(norm_expr))
filtered_expr_count <- norm_expr[keep_idx_expr, ]
cat( "Probes after filtering:", nrow(filtered_expr_count), "\n")
write.csv(as.data.frame(norm_expr), file = file.path(results_dir, "normalized_expression_all_probes.csv"))
write.csv(as.data.frame(filtered_expr_count), file = file.path(results_dir, "normalized_expression_filtered_probes.csv"))


## ---- Phenotype & groups ----
gse <- getGEO("GSE72267", GSEMatrix = TRUE)
gse_eset <- if (is.list(gse)) gse[[1]] else gse
pheno <- pData(gse_eset)
write.csv(pheno, file = file.path(results_dir, "phenotype_raw.csv"), row.names = TRUE)
# Heuristic column selection for group labels
found_col <- NULL
for (cname in c("source_name_ch1","title","characteristics_ch1","characteristics_ch1.1")) {
  if (cname %in% colnames(pheno)) { found_col <- cname; break }
}
if (is.null(found_col)) found_col <- colnames(pheno)[1]
raw_labels <- as.character(pheno[[found_col]])
group_label <- rep(NA, length(raw_labels))
group_label[grepl("parkinson|pd", raw_labels, ignore.case = TRUE)] <- "PD"
group_label[grepl("control|healthy|normal", raw_labels, ignore.case = TRUE)] <- "Control"
if (any(is.na(group_label)) && "characteristics_ch1" %in% colnames(pheno)) {
  ch <- as.character(pheno$characteristics_ch1)
  group_label[is.na(group_label) & grepl("Parkinson|PD", ch, ignore.case = TRUE)] <- "PD"
  group_label[is.na(group_label) & grepl("control|healthy|normal", ch, ignore.case = TRUE)] <- "Control"
}
group_label[is.na(group_label)] <- "Unknown"
groups <- factor(group_label, levels = c("Control","PD","Unknown"))
write.csv(data.frame(sample = rownames(pheno), group = groups), file = file.path(results_dir, "sample_groups.csv"), row.names = FALSE)
cat("Group counts:\n"); print(table(groups))

## ---- Differential expression with limma (PD vs Control) ----
# Remove Unknown samples
valid_idx <- which(groups %in% c("Control","PD"))
expr_for_de <- filtered_expr_count[, valid_idx]
group_for_de <- droplevels(groups[valid_idx])
design <- model.matrix(~0 + group_for_de)
colnames(design) <- levels(group_for_de)
fit <- lmFit(expr_for_de, design)
cont.matrix <- makeContrasts(PDvsControl = PD - Control, levels = design)
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2)
topTable_res <- topTable(fit2, coef = "PDvsControl", number = 200, adjust.method = "BH")
write.csv(topTable_res, file = file.path(results_dir, "limma_topTable_PDvsControl.csv"))

# Volcano plot
volcano_df <- data.frame(logFC = topTable_res$logFC, negLogP = -log10(topTable_res$P.Value), gene = rownames(topTable_res))
png(file.path(results_dir, "volcano_top200.png"), width = 900, height = 700)
plot(volcano_df$logFC, volcano_df$negLogP, pch = 20, main = "Volcano (top 200 genes)", xlab = "log2 FC", ylab = "-log10(P)")
dev.off()


## ---- Summary file ----
summary_file <- file.path(results_dir, "assignment_summary.txt")
cat("Assignment summary for dataset GSE72267\n", file = summary_file)
cat("Number of CEL files read: ", length(cel_files), "\n", file = summary_file, append = TRUE)
cat("Automatic raw outlier count (heuristic): ", raw_outlier_count_auto, "\n", file = summary_file, append = TRUE)
cat("Automatic normalized outlier count (heuristic): ", norm_outlier_count_auto, "\n", file = summary_file, append = TRUE)
cat("Probes before filtering:", nrow(norm_expr), "\n", file = summary_file, append = TRUE)
cat("Probes after filtering (expression-count method):", nrow(filtered_expr_count), "\n", file = summary_file, append = TRUE)
cat("Group counts (Control / PD / Unknown):\n", file = summary_file, append = TRUE)
capture.output(table(groups), file = summary_file, append = TRUE)
cat("arrayQualityMetrics (raw) dir: ", aqm_raw_dir, "\n", file = summary_file, append = TRUE)
cat("arrayQualityMetrics (normalized) dir: ", aqm_norm_dir, "\n", file = summary_file, append = TRUE)
cat("Top DE results (limma) saved to: ", file.path(results_dir, "limma_topTable_PDvsControl.csv"), "\n", file = summary_file, append = TRUE)
cat("All outputs written to:", results_dir, "\n", file = summary_file, append = TRUE)
cat("IMPORTANT: open the two index.html reports (raw & normalized) and visually report how many arrays are red-flagged for your assignment.\n", file = summary_file, append = TRUE)
