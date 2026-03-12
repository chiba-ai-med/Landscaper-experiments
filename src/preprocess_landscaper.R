source("src/Functions.R")

# Parameter
infile <- 'output/hpbase/integrated/seurat_annotated.RData'
outfile1 <- 'output/hpbase/integrated/seurat_annotated_landscaper.RData'
outfile2 <- 'output/hpbase/integrated/seurat.tsv'
outfile3 <- 'output/hpbase/integrated/group.tsv'
outfile4 <- 'output/hpbase/integrated_cov/group.tsv'
outfile5 <- 'output/hpbase/cont/group.tsv'
outfile6 <- 'output/hpbase/cont_cov/group.tsv'
outfile7 <- 'output/hpbase/DAPT/group.tsv'
outfile8 <- 'output/hpbase/DAPT_cov/group.tsv'
outfile9 <- 'output/hpbase/integrated/cov.tsv'
outfile10 <- 'output/hpbase/integrated_cov/cov.tsv'
outfile11 <- 'output/hpbase/cont/cov.tsv'
outfile12 <- 'output/hpbase/cont_cov/cov.tsv'
outfile13 <- 'output/hpbase/DAPT/cov.tsv'
outfile14 <- 'output/hpbase/DAPT_cov/cov.tsv'

# Loading
load(infile)

# Preprocess
## Only Ectoderm in 36h, 48h, 72h, 96h samples
target1 <- which(seurat.integrated@meta.data$germlayer == "Ectoderm")
target2 <- grep("36h|48h|72h|96h", seurat.integrated@meta.data$sample)
target <- intersect(target1, target2)
seurat.integrated <- seurat.integrated[, target]

## ---- ここを追加：sampleごとに2,751に揃える ----
set.seed(123456)

md <- seurat.integrated@meta.data
samp <- md$sample

idx_keep <- unlist(lapply(split(seq_len(ncol(seurat.integrated)), samp), function(ii) {
  if (length(ii) >= 2751) sample(ii, 2751) else ii
}), use.names = FALSE)

seurat.integrated <- seurat.integrated[, idx_keep]
## ----------------------------------------------

# Dimensional Reduction
## 純粋にNMFのみ（だめだったやり方）
# out <- as.matrix(seurat.integrated@assays$RNA@data[VariableFeatures(seurat.integrated), ])

# ## PCA（だめだったやり方）
# out <- seurat.integrated@reductions$pca@cell.embeddings[,1:20]

## PCA & UMAPやり直し
seurat.integrated <- RunPCA(seurat.integrated)
seurat.integrated <- RunUMAP(seurat.integrated, dims=1:30)
out <- seurat.integrated@reductions$pca@cell.embeddings[, 1:5]

# ランク5は分散の45％を説明
# > cumsum(seurat.integrated@reductions$pca@stdev^2) / sum(seurat.integrated@reductions$pca@stdev^2)
#  [1] 0.1811319 0.2855239 0.3489351 0.4059193 0.4519708 0.4917419 0.5285440
#  [8] 0.5629861 0.5901597 0.6151469 0.6361657 0.6558295 0.6745330 0.6927596
# [15] 0.7102934 0.7262952 0.7413836 0.7556817 0.7695486 0.7829664 0.7957339
# [22] 0.8076125 0.8184521 0.8290810 0.8394942 0.8494720 0.8591229 0.8682383
# [29] 0.8770848 0.8854015 0.8934235 0.9009822 0.9080725 0.9150929 0.9220168
# [36] 0.9284899 0.9346455 0.9407777 0.9464713 0.9519628 0.9573320 0.9625490
# [43] 0.9676759 0.9725660 0.9773622 0.9820812 0.9867006 0.9911952 0.9956385
# [50] 1.0000000

# set.seed(123456)
# out <- RunUMAP(seurat.integrated, dims=1:50, n.components=20L)@reductions$umap@cell.embeddings

# ## Normalization
# # Min-Max Normalization（だめだったやり方）
# # out <- out - min(out)
# out <- apply(out, 2, function(x){x - min(x)})

# Cluster Number
# group <- as.character(Idents(seurat.integrated))
group <- as.character(seurat.integrated@meta.data$celltype)

# Time
time <- unlist(.timevec[seurat.integrated@meta.data$sample])

# Output
save(seurat.integrated, file=outfile1)
write.table(out, outfile2, row.names=FALSE, col.names=FALSE, quote=FALSE)

write.table(group, outfile3, row.names=FALSE, col.names=FALSE, quote=FALSE)
write.table(group, outfile4, row.names=FALSE, col.names=FALSE, quote=FALSE)

write.table(group[grep("cont-", seurat.integrated@meta.data$sample)],
    outfile5, row.names=FALSE, col.names=FALSE, quote=FALSE)
write.table(group[grep("cont-", seurat.integrated@meta.data$sample)],
    outfile6, row.names=FALSE, col.names=FALSE, quote=FALSE)
write.table(group[grep("DAPT-", seurat.integrated@meta.data$sample)],
    outfile7, row.names=FALSE, col.names=FALSE, quote=FALSE)
write.table(group[grep("DAPT-", seurat.integrated@meta.data$sample)],
    outfile8, row.names=FALSE, col.names=FALSE, quote=FALSE)

write.table(time, outfile9, row.names=FALSE, col.names=FALSE, quote=FALSE)
write.table(time, outfile10, row.names=FALSE, col.names=FALSE, quote=FALSE)

write.table(time[grep("cont-", seurat.integrated@meta.data$sample)],
    outfile11, row.names=FALSE, col.names=FALSE, quote=FALSE)
write.table(time[grep("cont-", seurat.integrated@meta.data$sample)],
    outfile12, row.names=FALSE, col.names=FALSE, quote=FALSE)
write.table(time[grep("DAPT-", seurat.integrated@meta.data$sample)],
    outfile13, row.names=FALSE, col.names=FALSE, quote=FALSE)
write.table(time[grep("DAPT-", seurat.integrated@meta.data$sample)],
    outfile14, row.names=FALSE, col.names=FALSE, quote=FALSE)
