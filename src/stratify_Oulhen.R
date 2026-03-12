source("src/Functions.R")

# Parameter
infile1 <- commandArgs(trailingOnly=TRUE)[1]
infile2 <- commandArgs(trailingOnly=TRUE)[2]
outfile1 <- commandArgs(trailingOnly=TRUE)[3]
outfile2 <- commandArgs(trailingOnly=TRUE)[4]
outfile3 <- commandArgs(trailingOnly=TRUE)[5]
outfile4 <- commandArgs(trailingOnly=TRUE)[6]

# Loading
load(infile1)
seurat.integrated.original <- seurat.integrated
bindata <- read.table(infile2, header=FALSE)

# Stratification
idx_cont <- grep("Control", seurat.integrated.original@meta.data$orig.ident)
idx_dapt <- grep("DAPT", seurat.integrated.original@meta.data$orig.ident)

# Save（Seurat Object）
seurat.integrated <- seurat.integrated.original[, idx_cont]
save(seurat.integrated, file=outfile1)

seurat.integrated <- seurat.integrated.original[, idx_dapt]
save(seurat.integrated, file=outfile2)

# Save（BIN Data）
bindata_cont <- bindata[idx_cont, ]
bindata_dapt <- bindata[idx_dapt, ]

write.table(bindata_cont, outfile3, row.names=FALSE, col.names=FALSE, quote=FALSE)
write.table(bindata_dapt, outfile4, row.names=FALSE, col.names=FALSE, quote=FALSE)
