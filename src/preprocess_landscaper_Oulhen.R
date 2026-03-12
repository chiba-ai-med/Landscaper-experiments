source("src/Functions.R")

# Parameter
infile <- 'data/Oulhen/Combined_Dim20_res0.5.rds'
outfile1 <- 'output/oulhen/integrated/seurat_annotated_landscaper.RData'
outfile2 <- 'output/oulhen/integrated/seurat.tsv'
outfile3 <- 'output/oulhen/integrated/group.tsv'
outfile4 <- 'output/oulhen/cont/group.tsv'
outfile5 <- 'output/oulhen/DAPT/group.tsv'

# Loading
seurat.integrated <- readRDS(infile)

# Identsで2,4,8,10-12の細胞を除外
clusters_to_exclude <- c(2, 4, 8, 10, 11, 12)
cells_to_keep <- !Idents(seurat.integrated) %in% clusters_to_exclude
seurat.integrated <- seurat.integrated[, cells_to_keep]

# Cell-type Annotation（クラスタ番号を細胞型名に変換）
celltype_mapping <- c(
  "0" = "Aboral_ectoderm",
  "1" = "Ciliary_band",
  "3" = "Oral_ectoderm",
  "5" = "Aboral_ectoderm",
  "6" = "Ciliary_band",
  "7" = "Ciliary_band",
  "9" = "Neurons"
)

# 新しいメタデータ列を追加
seurat.integrated@meta.data$celltype <- celltype_mapping[as.character(Idents(seurat.integrated))]

# 色の定義（RGBから16進数に変換）
celltype_colors_map <- c(
  "Aboral_ectoderm" = "#008080",  # RGB(0, 128, 128) - ティール
  "Oral_ectoderm"   = "#FFFF00",  # RGB(255, 255, 0) - イエロー
  "Neurons"         = "#FF00FF",  # RGB(255, 0, 255) - マゼンタ
  "Ciliary_band"     = "#00FF26"   # RGB(0, 255, 38) - グリーン
)

# 各細胞に対応する色を割り当て
seurat.integrated@meta.data$celltype_colors <- celltype_colors_map[seurat.integrated@meta.data$celltype]

# celltypeをIdentsに設定（オプション）
Idents(seurat.integrated) <- seurat.integrated@meta.data$celltype

# Dimensional Reduction
# 古いPCAを削除して再計算
seurat.integrated@reductions$pca <- NULL  # 既存のPCAを削除

# Variable featuresを再計算（必要な場合）
seurat.integrated <- FindVariableFeatures(seurat.integrated, selection.method = "vst", nfeatures = 2000)

# データをスケーリング
seurat.integrated <- ScaleData(seurat.integrated)

# PCAを新しく計算
seurat.integrated <- RunPCA(seurat.integrated, features = VariableFeatures(seurat.integrated))

# 確認
print(paste("New PCA dimensions:", nrow(seurat.integrated@reductions$pca@cell.embeddings), "cells x", 
            ncol(seurat.integrated@reductions$pca@cell.embeddings), "PCs"))

# UMAPも再計算
seurat.integrated@reductions$umap <- NULL  # 既存のUMAPも削除
set.seed(123456)
seurat.integrated <- RunUMAP(seurat.integrated, dims=1:30)
# 出力データの準備
out <- seurat.integrated@reductions$pca@cell.embeddings[, 1:12]

# ランク6は分散の39％を説明
# > cumsum(seurat.integrated@reductions$pca@stdev^2) / sum(seurat.integrated@reductions$pca@stdev^2)
# [1] 0.1716118 0.2311870 0.2808243 0.3296187 0.3645635 0.3990872 0.4252697
#  [8] 0.4458435 0.4654341 0.4842250 0.5017786 0.5183489 0.5340991 0.5494892
# [15] 0.5647924 0.5793532 0.5934815 0.6074565 0.6212074 0.6348759 0.6482460
# [22] 0.6615662 0.6746894 0.6876892 0.7006632 0.7134594 0.7261019 0.7386621
# [29] 0.7511957 0.7636043 0.7758493 0.7880699 0.8002443 0.8123642 0.8244212
# [36] 0.8364300 0.8483877 0.8602890 0.8721618 0.8839772 0.8957385 0.9074491
# [43] 0.9191520 0.9307913 0.9424104 0.9539901 0.9655381 0.9770580 0.9885455
# [50] 1.0000000

# Cluster Number
group <- as.character(seurat.integrated@meta.data$celltype)

# Output
save(seurat.integrated, file=outfile1)
write.table(out, outfile2, row.names=FALSE, col.names=FALSE, quote=FALSE)

write.table(group, outfile3, row.names=FALSE, col.names=FALSE, quote=FALSE)
write.table(group[grep("Control", seurat.integrated@meta.data$orig.ident)],
    outfile4, row.names=FALSE, col.names=FALSE, quote=FALSE)
write.table(group[grep("DAPT", seurat.integrated@meta.data$orig.ident)],
    outfile5, row.names=FALSE, col.names=FALSE, quote=FALSE)
