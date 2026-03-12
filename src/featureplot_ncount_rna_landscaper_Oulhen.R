source("src/Functions.R")

# Parameter
infile <- commandArgs(trailingOnly=TRUE)[1]
outfile <- commandArgs(trailingOnly=TRUE)[2]

# Loading
load(infile)

# Log10
seurat.integrated[["Log10_nCount_RNA"]] <- log10(seurat.integrated@meta.data$nCount_RNA)

# SeuratからSingleCellExperimentに変換
seurat.integrated <- JoinLayers(seurat.integrated)
sce <- as.SingleCellExperiment(seurat.integrated)

# colDataのカラム名を確認（デバッグ用、必要に応じてコメントアウト）
# print(colnames(colData(sce)))

# Hexbinの計算
sce <- make_hexbin(sce, 
                   nbins = 28,
                   dimension_reduction = "UMAP")

# Plot
png(file=outfile, width=600, height=600)
plot_hexbin_meta(sce, 
                 col = "Log10_nCount_RNA",
                 action = "mean") +
    xlim(c(-15,15)) + ylim(c(-15,15)) + 
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