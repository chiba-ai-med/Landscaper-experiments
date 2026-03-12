source("src/Functions.R")

# Parameter
infile <- commandArgs(trailingOnly=TRUE)[1]
outfile <- commandArgs(trailingOnly=TRUE)[2]

# Loading
load(infile)

# sample名から時間情報を抽出して新しい列を作成
seurat.integrated$timepoint <- gsub("cont-|DAPT-", "", seurat.integrated$sample)

# 見やすい赤系、青系、緑系の色を定義
time_colors <- c(
  "24h" = "#E74C3C",  # 赤系
  "36h" = "#3498DB",  # 青系
  "48h" = "#27AE60",   # 緑系
  "72h" = "#F39C12",  # オレンジ系
  "96h" = "#8E44AD"   # 紫系
)

# Plot
png(file=outfile, width=600, height=600)
g <- DimPlot(seurat.integrated, 
             reduction = "umap", 
             group.by = "timepoint", 
             cols = time_colors,
             label = TRUE, 
             pt.size = 2, 
             label.size = 0) + 
    NoLegend() +
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
print(g)
dev.off()