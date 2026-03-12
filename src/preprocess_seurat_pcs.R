source("src/Functions.R")

# Arguments
args <- commandArgs(trailingOnly=TRUE)
infile <- args[1]
outfile <- args[2]
pc <- as.integer(args[3])

# Loading
load(infile)

# Preprocess
out <- seurat.integrated@reductions$pca@cell.embeddings[, 1:pc]

# Output
write.table(out, outfile, row.names=FALSE, col.names=FALSE, quote=FALSE)
