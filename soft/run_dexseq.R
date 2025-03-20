# Load required libraries
library(DEXSeq)
library(BiocParallel)

# Function to parse group file
parse_group_file <- function(file_path) {
  lines <- readLines(file_path)
  sample_files <- c()
  conditions <- c()
  for (line in lines) {
    parts <- strsplit(line, ":")[[1]]
    group_name <- parts[1]
    files <- strsplit(parts[2], ",")[[1]]
    sample_files <- c(sample_files, files)
    conditions <- c(conditions, rep(group_name, length(files)))
  }
  list(sample_files = sample_files, conditions = conditions)
}

# Function to safely convert DEXSeq results to a data frame
safe_dexseq_to_df <- function(dxr) {
  # Convert to data frame
  df <- as.data.frame(dxr)
  
  # Check for list columns and convert them
  for (col in names(df)) {
    if (is.list(df[[col]])) {
      # Convert list column to character strings
      df[[col]] <- sapply(df[[col]], function(x) {
        if (length(x) > 1) {
          paste(x, collapse = ";")
        } else if (length(x) == 1) {
          as.character(x)
        } else {
          NA
        }
      })
    }
  }
  
  return(df)
}

# Get command-line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  stop("Error: Please provide the paths to: 1) group file, 2) DEXSeq-prepared annotation file, 3) output directory")
}

group_file <- args[1]       # First argument is the group file path
dexseq_gff_file <- args[2]  # Second argument is the DEXSeq-prepared annotation file
output_dir <- args[3]       # Third argument is the output directory

# Create output directory if it doesn't exist
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Check if DEXSeq-prepared annotation file exists
if (!file.exists(dexseq_gff_file)) {
  message("DEXSeq-prepared annotation file not found.")
  message("Please run dexseq_prepare_annotation.py first:")
  message("python /path/to/dexseq_prepare_annotation.py your_gtf_file.gtf ", dexseq_gff_file)
  stop("Annotation file missing.")
}

# Parse the group file
sample_info <- parse_group_file(group_file)
sample_files <- sample_info$sample_files
conditions <- sample_info$conditions

# Create sample table
sampleTable <- data.frame(
  row.names = basename(sample_files),
  condition = factor(conditions),
  stringsAsFactors = FALSE
)

# Set up parallel processing (optional but recommended for speed)
BPPARAM <- MulticoreParam(workers = 4)

# Step 1: Create DEXSeqDataSet from count files
message("Creating DEXSeqDataSet from count files...")
dxd <- DEXSeq::DEXSeqDataSetFromHTSeq(
  countfiles = sample_files,
  sampleData = sampleTable,
  design = ~ sample + exon + condition:exon,
  flattenedfile = dexseq_gff_file
)

# Step 2: Normalization
message("Performing normalization...")
dxd <- DEXSeq::estimateSizeFactors(dxd)

# Step 3: Dispersion estimation
message("Estimating dispersions...")
dxd <- DEXSeq::estimateDispersions(dxd, BPPARAM = BPPARAM)

# Step 4: Testing for differential exon usage
message("Testing for differential exon usage...")
dxd <- DEXSeq::testForDEU(dxd, BPPARAM = BPPARAM)

# Step 5: Estimate fold changes
message("Estimating fold changes...")
# Use the first condition as reference level
dxd$condition <- relevel(dxd$condition, ref = levels(dxd$condition)[1])
dxd <- DEXSeq::estimateExonFoldChanges(dxd, fitExpToVar = "condition")

# Step 6: Extract results
message("Extracting results...")
dxr <- DEXSeq::DEXSeqResults(dxd)

# Step 7: Save results to files
message("Saving results...")
full_results_file <- file.path(output_dir, "DEXSeq_full_results.csv")
significant_results_file <- file.path(output_dir, "DEXSeq_significant_results.csv")
summary_file <- file.path(output_dir, "DEXSeq_summary.csv")

# Save full results table - using our safe conversion function
results_df <- safe_dexseq_to_df(dxr)
write.csv(results_df, file = full_results_file, row.names = FALSE)

# Save significant results (adjusted p-value < 0.05)
significant_idx <- which(dxr$padj < 0.05)
if (length(significant_idx) > 0) {
  dxr_significant <- dxr[significant_idx,]
  significant_df <- safe_dexseq_to_df(dxr_significant)
  write.csv(significant_df, file = significant_results_file, row.names = FALSE)
} else {
  # Create an empty file if no significant results
  message("No significant results found at adjusted p-value < 0.05")
  file.create(significant_results_file)
}

# Create and save a summary of results
# Extract the condition name for the fold change column dynamically
condition_levels <- levels(dxd$condition)
fold_change_col <- paste0("log2fold_", condition_levels[2], "_", condition_levels[1])

# Make sure the fold change column exists in our results
if (fold_change_col %in% colnames(results_df)) {
  up_regulated <- sum(results_df$padj < 0.05 & results_df[[fold_change_col]] > 0, na.rm = TRUE)
  down_regulated <- sum(results_df$padj < 0.05 & results_df[[fold_change_col]] < 0, na.rm = TRUE)
} else {
  # If column doesn't exist, set counts to NA
  message("Fold change column not found in results")
  up_regulated <- NA
  down_regulated <- NA
}

summary <- data.frame(
  Total_exons = nrow(results_df),
  Significant_exons = sum(results_df$padj < 0.05, na.rm = TRUE),
  Up_regulated = up_regulated,
  Down_regulated = down_regulated,
  Reference_condition = condition_levels[1],
  Test_condition = condition_levels[2]
)
write.csv(summary, file = summary_file, row.names = FALSE)

# Optional: Create a basic plot for visualization
ma_plot_file <- file.path(output_dir, "DEXSeq_MA_plot.pdf")
pdf(ma_plot_file, width = 10, height = 8)
plotMA(dxr, cex = 0.8)
dev.off()

# Optional: Plot top significant genes
if (sum(results_df$padj < 0.05, na.rm = TRUE) > 0) {
  # Get top 5 genes (or fewer if less are available)
  sig_genes <- unique(results_df$groupID[results_df$padj < 0.05 & !is.na(results_df$padj)])
  top_genes <- sig_genes[1:min(5, length(sig_genes))]
  
  for (gene in top_genes) {
    gene_plot_file <- file.path(output_dir, paste0("DEXSeq_", gene, "_plot.pdf"))
    pdf(gene_plot_file, width = 12, height = 8)
    plotDEXSeq(dxr, gene, legend = TRUE, cex.axis = 1.2, cex = 1.3, lwd = 2)
    dev.off()
  }
}

# Session info for reproducibility
session_info_file <- file.path(output_dir, "DEXSeq_session_info.txt")
capture.output(sessionInfo(), file = session_info_file)

message("Analysis complete! Results saved to ", output_dir)