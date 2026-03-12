source("src/Functions.R")

# Parameter
infile1 <- commandArgs(trailingOnly=TRUE)[1]
outfile1 <- commandArgs(trailingOnly=TRUE)[2]
outfile2 <- commandArgs(trailingOnly=TRUE)[3]

# Loading
load(infile1)
seurat.integrated.original <- seurat.integrated

# Stratification
idx_cont <- grep("Control", seurat.integrated.original@meta.data$orig.ident)
idx_dapt <- grep("DAPT", seurat.integrated.original@meta.data$orig.ident)

# Save (Seurat Object)
seurat.integrated <- seurat.integrated.original[, idx_cont]
save(seurat.integrated, file=outfile1)

seurat.integrated <- seurat.integrated.original[, idx_dapt]
save(seurat.integrated, file=outfile2)
