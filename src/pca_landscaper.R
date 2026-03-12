source("src/Functions.R")

# Parameter
args <- commandArgs(trailingOnly=TRUE)
infile <- args[1]  # seurat_annotated_landscaper.RData
outfile <- args[2]  # pca.png
# infile <- 'output/hpbase/integrated/seurat_annotated_landscaper.RData'
# outfile <- 'plot/hpbase/integrated/pca.png'

# Load
load(infile)

# PCA Pair Plot
coordinates <- as.data.frame(seurat.integrated@reductions$pca@cell.embeddings)
celltypes <- seurat.integrated@meta.data$celltype

# 色の定義（RGBから16進数に変換）
celltype_colors_map <- c(
  "Aboral_ectoderm" = "#008080",  # RGB(0, 128, 128) - ティール
  "Oral_ectoderm"   = "#FFFF00",  # RGB(255, 255, 0) - イエロー
  "Neurons"         = "#FF00FF",  # RGB(255, 0, 255) - マゼンタ
  "Ciliary_band"     = "#00FF26"   # RGB(0, 255, 38) - グリーン
)

# Plot
png(file = outfile, width = 1200, height = 1200, res = 150)
pairs(
  coordinates[, 1:15],
  col = celltype_colors_map[celltypes],
  pch = 16,
  cex = 0.5,
  main = "PCA Pair Plot"
)
dev.off()
