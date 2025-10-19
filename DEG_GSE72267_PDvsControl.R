# =====================================================================
#             AI and Biotechnology / Bioinformatics Internship 2025
# =====================================================================
#       Module II: Microarray Data Analysis — Assignment (GSE72267)
# =====================================================================
#                Parkinson’s Disease vs Control Samples
# =====================================================================

gc()  # Clear memory

#### 1. Install and Load Packages ####
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install(c("limma", "AnnotationDbi", "hgu133a2.db"), ask = FALSE)
install.packages(c("dplyr", "tibble", "ggplot2", "pheatmap"), dependencies = TRUE)

library(limma)
library(AnnotationDbi)
library(hgu133a2.db)    # For HG-U133A_2 platform
library(dplyr)
library(tibble)
library(ggplot2)
library(pheatmap)

#### 2. Load Processed Data ####
# These files should come from your preprocessing script
results_dir <- file.path(getwd(), "Results")
expr_file   <- file.path(results_dir, "normalized_expression_filtered_probes.csv")
group_file  <- file.path(results_dir, "sample_groups.csv")

expr_df <- read.csv(expr_file, row.names = 1, check.names = FALSE)
sample_groups <- read.csv(group_file, stringsAsFactors = FALSE)

expr_mat <- as.matrix(expr_df)

# Remove suffixes like "_001MBOC.CEL.gz" and keep only GSM IDs
# Reload sample_groups from CSV
sample_groups <- read.csv("Results/sample_groups.csv", stringsAsFactors = FALSE)

# Check its structure
head(sample_groups)

# Set rownames properly
rownames(sample_groups) <- sample_groups$sample

# Clean the column names of the expression matrix
colnames(expr_mat) <- sub("_.*", "", colnames(expr_mat))

# Find common samples
common_samples <- intersect(colnames(expr_mat), sample_groups$sample)
cat("Common samples found:", length(common_samples), "\n")

# Subset both objects
expr_mat <- expr_mat[, common_samples, drop = FALSE]
sample_groups <- sample_groups[common_samples, , drop = FALSE]

# Confirm alignment
cat("Expression matrix dimensions:", dim(expr_mat), "\n")
head(sample_groups)


#### 3. Probe → Gene Mapping ####
probe_ids <- rownames(expr_mat)

gene_symbols <- mapIds(
  hgu133a2.db,
  keys = probe_ids,
  keytype = "PROBEID",
  column = "SYMBOL",
  multiVals = "first"
)

mapping_df <- data.frame(PROBEID = probe_ids, SYMBOL = gene_symbols, stringsAsFactors = FALSE)
mapping_df <- mapping_df %>% filter(!is.na(SYMBOL))

expr_mapped <- expr_mat[mapping_df$PROBEID, , drop = FALSE]
id_vector <- mapping_df$SYMBOL

# Collapse duplicate probes (average by gene)
expr_gene <- limma::avereps(expr_mapped, ID = id_vector)

cat("Genes after collapsing duplicates:", nrow(expr_gene), "\n")

#### 4. Define Groups ####
# Normalize group names
sample_groups$group <- ifelse(grepl("control", sample_groups$group, ignore.case = TRUE), "Control",
                              ifelse(grepl("pd|parkinson", sample_groups$group, ignore.case = TRUE), "PD",
                                     sample_groups$group))

groups <- factor(sample_groups$group, levels = c("Control", "PD"))
table(groups)

# Keep only PD and Control
valid_idx <- which(groups %in% c("Control", "PD"))
expr_de <- expr_gene[, valid_idx, drop = FALSE]
grp_de <- droplevels(groups[valid_idx])

#### 5. Linear Modeling (limma) ####
design <- model.matrix(~0 + grp_de)
colnames(design) <- levels(grp_de)

fit <- lmFit(expr_de, design)
contrast_matrix <- makeContrasts(PD_vs_Control = PD - Control, levels = design)
fit2 <- contrasts.fit(fit, contrast_matrix)
fit2 <- eBayes(fit2)

#### 6. Differential Expression ####
deg_results <- topTable(fit2, coef = "PD_vs_Control", number = Inf, adjust.method = "BH")

deg_results$Regulation <- ifelse(
  deg_results$adj.P.Val < 0.05 & deg_results$logFC > 1, "Upregulated",
  ifelse(deg_results$adj.P.Val < 0.05 & deg_results$logFC < -1, "Downregulated", "Not Significant")
)

upregulated <- subset(deg_results, Regulation == "Upregulated")
downregulated <- subset(deg_results, Regulation == "Downregulated")

cat("Upregulated:", nrow(upregulated), " | Downregulated:", nrow(downregulated), "\n")

dir.create("Results_DEG", showWarnings = FALSE)
write.csv(deg_results, "Results_DEG/All_DEGs_PD_vs_Control.csv")
write.csv(upregulated, "Results_DEG/Upregulated_DEGs.csv")
write.csv(downregulated, "Results_DEG/Downregulated_DEGs.csv")

#### 7. Visualization ####
dir.create("Result_Plots", showWarnings = FALSE)

# --- Volcano Plot ---
png("Result_Plots/Volcano_PD_vs_Control.png", width = 2000, height = 1500, res = 300)
ggplot(deg_results, aes(x = logFC, y = -log10(adj.P.Val), color = Regulation)) +
  geom_point(alpha = 0.7, size = 2) +
  scale_color_manual(values = c("Upregulated" = "red",
                                "Downregulated" = "blue",
                                "Not Significant" = "grey")) +
  theme_minimal() +
  labs(title = "Volcano Plot: Parkinson’s Disease vs Control",
       x = "log2 Fold Change", y = "-log10(Adjusted P-value)",
       color = "Regulation")
dev.off()

# --- Heatmap of Top 25 DEGs ---
top_genes <- head(rownames(deg_results[order(deg_results$adj.P.Val), ]), 25)
heatmap_data <- expr_de[top_genes, ]

png("Result_Plots/Heatmap_Top25_PDvsControl.png", width = 2000, height = 1500, res = 300)
pheatmap(heatmap_data,
         scale = "row",
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         show_rownames = TRUE,
         show_colnames = TRUE,
         main = "Top 25 DEGs (PD vs Control)",
         color = colorRampPalette(c("blue", "white", "red"))(100))
dev.off()

#### 8. Short Summary ####
cat("
Result Summary:
- Dataset: GSE72267 (Parkinson’s Disease vs Control)
- Platform: Affymetrix HG-U133A_2
- Probes mapped to gene symbols using hgu133a2.db
- Duplicate probes collapsed by averaging (limma::avereps)
- Contrast: PD_vs_Control
- Significant DEGs: adj.P.Val < 0.05 and |log2FC| > 1
- Upregulated genes:", nrow(upregulated), " Downregulated genes:", nrow(downregulated), "
Results saved to: Results_DEG and Result_Plots folders.
")

