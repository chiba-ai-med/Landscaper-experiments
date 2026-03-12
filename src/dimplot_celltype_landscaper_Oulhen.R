source("src/Functions.R")

# Parameter
infile <- commandArgs(trailingOnly=TRUE)[1]
outfile <- commandArgs(trailingOnly=TRUE)[2]

# Loading
load(infile)

# 色
cols <- sapply(names(table(seurat.integrated$celltype)), function(x){
    target <- which(seurat.integrated$celltype == x)[1]
    tmp <- seurat.integrated$celltype_colors[target]
    names(tmp) <- seurat.integrated$celltype[x]
    tmp
})
names(cols) <- gsub("\\.NA", "", names(cols))

# Idents切り替え
Idents(seurat.integrated) <- factor(seurat.integrated$celltype, levels=names(cols))

# Plot
png(file=outfile, width=600, height=600)
DimPlot(seurat.integrated, reduction = "umap", label=FALSE, pt.size=3, label.size=6, cols=cols) + NoLegend() +
    theme(axis.line = element_blank(),
           axis.text.x = element_blank(),
           axis.text.y = element_blank(),
           axis.ticks = element_blank(),
           axis.title.x = element_blank(),
           axis.title.y = element_blank(),
           panel.background = element_blank(),
           panel.border = element_blank(),
           panel.grid.major = element_blank(),
           panel.grid.minor = element_blank(),
           plot.background = element_blank())
dev.off()