#----------------------------------------------------------------------------------------------------
####          Assignment 2 ####
# ----------------------------------------------------------------------------------------------
# Function to classify genes
classify_gene <- function(logFC, padj) {
  if (is.na(padj)) padj <- 1  # Replace missing padj with 1
  if (logFC > 1 && padj < 0.05) {
    return("Upregulated")
  } else if (logFC < -1 && padj < 0.05) {
    return("Downregulated")
  } else {
    return("Not_Significant")
  }
}
#--------------------------------------------------------------------------------------------------
# Directories and files

input_dir <- "C:/Users/PMLS/Desktop/AI_Omics_Internship_2025/Raw_data"    # Replace with actual directory
output_dir <- file.path(input_dir, "Results")
if (!dir.exists(output_dir)) dir.create(output_dir)

files_to_process <- c("DEGs_data_1.csv", "DEGs_data_2.csv")
# Loop through files
for (file_name in files_to_process) {

  cat("\nProcessing:", file_name, "\n")
  
  # Read file
  input_file_path <- file.path(input_dir, file_name)
  data <- read.csv(input_file_path, stringsAsFactors = FALSE)
  
  # Replace missing padj with 1
  data$padj[is.na(data$padj)] <- 1
  
  # Apply classify_gene() function to create 'status' column
  data$status <- mapply(classify_gene, data$logFC, data$padj)
  # Save processed file
  output_file_path <- file.path(output_dir, paste0("processed_", file_name))
  write.csv(data, output_file_path, row.names = FALSE)
  cat("Processed file saved to:", output_file_path, "\n")
  # Print summary counts
  cat("\nSummary for", file_name, ":\n")
  print(table(data$status))
  
}
save.image("MaryamAhmedAli_Class_2_Assignment.RData")

# -----------------------------------------------------------------------------------------