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

# Run DESeq2
dds <- DESeqDataSetFromMatrix(countData = round(count_mat),
                              colData = metadata,
                              design = ~ condition)
dds <- DESeq(dds)
res <- results(dds)

# Save results
saveRDS(dds, snakemake@output[["dds_rds"]])
write.csv(as.data.frame(res), snakemake@output[["results_csv"]])
