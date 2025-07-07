library(DESeq2)

# Read metadata
metadata <- read.csv(snakemake@input[["metadata"]], row.names = 1)

# Match sample names to abundance files
files <- snakemake@input[["abundance"]]
names(files) <- rownames(metadata)

# Read abundance files and extract estimated counts using base R
count_list <- lapply(files, function(f) {
  tab <- read.table(f, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  setNames(tab$est_counts, tab$target_id)
})

# Merge into a single count matrix
all_genes <- unique(unlist(lapply(count_list, names)))
count_mat <- sapply(count_list, function(x) {
  x <- x[match(all_genes, names(x))]
  x[is.na(x)] <- 0
  x
})
rownames(count_mat) <- all_genes

# CRITICAL: Check and set reference level explicitly
print("Available conditions:")
print(unique(metadata$condition))

print("Current factor levels:")
print(levels(as.factor(metadata$condition)))

# DESeq2 uses ALPHABETICAL order by default for reference level
# The FIRST level alphabetically becomes the reference (control)
print("DESeq2 will use this as reference (control) by default:")
print(sort(unique(metadata$condition))[1])

# Let DESeq2 use alphabetical order (current behavior)
# Uncomment these lines if you want to keep current behavior:
dds <- DESeqDataSetFromMatrix(countData = round(count_mat),
                              colData = metadata,
                              design = ~ condition)

# OPTION 2: Explicitly set your desired reference level
# Replace "your_control_condition" with your actual control group name
# For example, if your conditions are "treated" and "control", use "control"
# metadata$condition <- factor(metadata$condition, 
                            # levels = c("control", "treated"))  # Put control first

# Alternative way to set reference if you know the control name:
# metadata$condition <- relevel(factor(metadata$condition), ref = "control")

# dds <- DESeqDataSetFromMatrix(countData = round(count_mat),
#                               colData = metadata,
#                               design = ~ condition)

# Verify the reference level
print("Final reference level (control):")
print(levels(dds$condition)[1])

# Run DESeq2
dds <- DESeq(dds)

# Get results - this compares all other levels to the reference
res <- results(dds)

# You can also get specific contrasts:
# res <- results(dds, contrast = c("condition", "treated", "control"))

print("Results comparison:")
print(paste("Positive log2FoldChange = higher in", 
            levels(dds$condition)[2], "vs", levels(dds$condition)[1]))
print(paste("Negative log2FoldChange = higher in", 
            levels(dds$condition)[1], "vs", levels(dds$condition)[2]))

# Save results
saveRDS(dds, snakemake@output[["dds_rds"]])
write.csv(as.data.frame(res), snakemake@output[["results_csv"]])

# Save a summary file with interpretation
interpretation <- data.frame(
  reference_level = levels(dds$condition)[1],
  comparison_level = levels(dds$condition)[2],
  positive_FC_means = paste("Higher in", levels(dds$condition)[2]),
  negative_FC_means = paste("Higher in", levels(dds$condition)[1]),
  total_genes = nrow(res),
  significant_genes = sum(res$padj < 0.05, na.rm = TRUE)
)

write.csv(interpretation, snakemake@output[["interp_csv"]], row.names = FALSE)