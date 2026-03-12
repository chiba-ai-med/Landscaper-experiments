source("src/Functions.R")

# Parameter
pc <- commandArgs(trailingOnly=TRUE)[1]
outfile1 <- commandArgs(trailingOnly=TRUE)[2]
outfile2 <- commandArgs(trailingOnly=TRUE)[3]
outfile3 <- commandArgs(trailingOnly=TRUE)[4]
outfile4 <- commandArgs(trailingOnly=TRUE)[5]
outfile5 <- commandArgs(trailingOnly=TRUE)[6]
outfile6 <- commandArgs(trailingOnly=TRUE)[7]

# Loading
all_states_DAPT <- unlist(read.delim(paste0('plot/oulhen/DAPT/', pc, '/Landscaper/Allstates.tsv'), header=FALSE, sep="|"))
bin_data_DAPT <- unlist(read.delim(paste0('output/oulhen/DAPT/', pc, '/binpca/BIN_DATA.tsv'), header=FALSE, sep="|"))
energy_DAPT <- unlist(read.table(paste0('plot/oulhen/DAPT/', pc, '/Landscaper/E.tsv'), header=FALSE))
load('output/oulhen/DAPT/seurat_annotated_landscaper.RData')

# Sort
names(all_states_DAPT) <- NULL
names(bin_data_DAPT) <- NULL

bin_data_DAPT <- gsub("\t", " ", bin_data_DAPT)
target_DAPT <- match(bin_data_DAPT, all_states_DAPT)

energy_DAPT_sorted <- energy_DAPT[target_DAPT]
energy_DAPT_sorted <- as.matrix(energy_DAPT_sorted)
rownames(energy_DAPT_sorted) <- colnames(seurat.integrated)

# Assign Labels
seurat.integrated$energy <- - energy_DAPT_sorted

# Seurat => SingleCellExperiment
sce <- SingleCellExperiment(assays = list(counts = seurat.integrated@assays$RNA$scale.data))
reducedDims(sce) <- list(UMAP = seurat.integrated@reductions$umap@cell.embeddings)
logcounts(sce) <- log10(counts(sce) - min(counts(sce)) + 1)

# Assign Energy
colData(sce)$energy <- seurat.integrated$energy

# Setting Hex bin
sce <- make_hexbin(sce, nbins = 50, dimension_reduction = "UMAP")
save(sce, file=outfile1)

# Plot
g <- FeaturePlot(seurat.integrated, reduction = "umap", features="energy", pt.size=1) + scale_colour_gradientn(colours = viridis(100)) + ggtitle("- Energy")
ggsave(file=outfile2, g, dpi=200, width=6, height=6)

# Plot
g <- plot_hexbin_meta(sce, col = "energy", action = "median", title = "") + ggtitle("- Energy")
ggsave(file=outfile3, g, dpi=200, width=6, height=6)

# Plot
g <- FeaturePlot(seurat.integrated, reduction = "umap", features="energy", pt.size=1) + scale_colour_gradientn(colours = viridis(100), limits=energy_limits) + ggtitle("- Energy")
ggsave(file=outfile4, g, dpi=200, width=6, height=6)

# Plot
g <- plot_hexbin_meta(sce, col = "energy", action = "median", title = "") + 
     scale_fill_gradientn(colours = viridis(100), 
                          limits = energy_limits,
                          oob = scales::squish) +  # 範囲外の値を圧縮
     ggtitle("- Energy")
ggsave(file=outfile5, g, dpi=200, width=6, height=6)

# 新規追加: エネルギーのマイナス値の等高線プロット（黄色ベース）
df <- as.data.frame(seurat.integrated@reductions$umap@cell.embeddings)
colnames(df) <- c("UMAP_1", "UMAP_2")
df$neg_energy <- seurat.integrated$energy  # エネルギーのマイナス値

g <- ggplot(df, aes(x = UMAP_1, y = UMAP_2)) +
  geom_point(color = "grey60", size = 0.5, alpha = 0.4) +
  stat_density_2d(aes(fill = ..density.., z = neg_energy), 
                  geom = "tile", 
                  contour = FALSE,
                  n = 100,
                  alpha = 0.45) +
  scale_fill_gradient(low = "white", high = "#FFB300") +  # 黄色系のグラデーション
  geom_density_2d(aes(z = neg_energy), 
                  color = "#CC8800", 
                  alpha = 0.7, 
                  size = 0.4) +  # 濃い黄色の等高線
  theme_void() +
  theme(legend.position = 'none',
        plot.background = element_rect(fill = "white", color = NA))

ggsave(file=outfile6, g, dpi=200, width=6, height=6)