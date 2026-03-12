source("src/Functions.R")

# Parameter
infile1 <- commandArgs(trailingOnly=TRUE)[1]  # integrated seurat_annotated_landscaper.RData
infile2 <- commandArgs(trailingOnly=TRUE)[2]  # integrated BIN_DATA.tsv
outfile1 <- commandArgs(trailingOnly=TRUE)[3] # cont BIN_DATA.tsv
outfile2 <- commandArgs(trailingOnly=TRUE)[4] # DAPT BIN_DATA.tsv

# Loading
load(infile1)
bindata <- read.table(infile2, header=FALSE)

# Stratification
idx_cont <- grep("Control", seurat.integrated@meta.data$orig.ident)
idx_dapt <- grep("DAPT", seurat.integrated@meta.data$orig.ident)

# Save (BIN Data)
write.table(bindata[idx_cont, ], outfile1, row.names=FALSE, col.names=FALSE, quote=FALSE)
write.table(bindata[idx_dapt, ], outfile2, row.names=FALSE, col.names=FALSE, quote=FALSE)
