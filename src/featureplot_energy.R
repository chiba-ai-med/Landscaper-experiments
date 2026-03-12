source("src/Functions.R")

args <- commandArgs(trailingOnly = TRUE)
input_rdata <- args[1]
outfile1 <- args[2]  # sce.RData
outfile2 <- args[3]  # energy.png
outfile3 <- args[4]  # energy_hex.png
outfile4 <- args[5]  # energy_rescaled.png
outfile5 <- args[6]  # energy_rescaled_hex.png
outfile6 <- args[7]  # energy_splitby.png
outfile7 <- args[8]  # energy_contour.png

# Load enriched Seurat
load(input_rdata)

# Negate energy for plotting (- Energy)
seurat.integrated$energy <- -seurat.integrated$energy

# Seurat => SingleCellExperiment
sce <- SingleCellExperiment(assays = list(counts = seurat.integrated@assays$RNA@counts))
reducedDims(sce) <- list(UMAP = seurat.integrated@reductions$umap@cell.embeddings)
logcounts(sce) <- log10(counts(sce) + 1)

# Assign Energy
colData(sce)$energy <- seurat.integrated$energy

# Setting Hex bin
sce <- make_hexbin(sce, nbins = 50, dimension_reduction = "UMAP")
save(sce, file=outfile1)

# Plot: energy
g <- FeaturePlot(seurat.integrated, reduction = "umap", features="energy", pt.size=1) + scale_colour_gradientn(colours = viridis(100)) + ggtitle("- Energy")
ggsave(file=outfile2, g, dpi=200, width=6, height=6)

# Plot: energy hex
g <- plot_hexbin_meta(sce, col = "energy", action = "median", title = "") + ggtitle("- Energy")
ggsave(file=outfile3, g, dpi=200, width=6, height=6)

# Plot: energy rescaled
g <- FeaturePlot(seurat.integrated, reduction = "umap", features="energy", pt.size=1) + scale_colour_gradientn(colours = viridis(100), limits=energy_limits) + ggtitle("- Energy")
ggsave(file=outfile4, g, dpi=200, width=6, height=6)

# Plot: energy rescaled hex
g <- plot_hexbin_meta(sce, col = "energy", action = "median", title = "") +
     scale_fill_gradientn(colours = viridis(100),
                          limits = energy_limits,
                          oob = scales::squish) +
     ggtitle("- Energy")
ggsave(file=outfile5, g, dpi=200, width=6, height=6)

# Plot: energy splitby
n_samples <- length(unique(seurat.integrated$sample))
plot_width <- max(750 * n_samples, 3000)
g <- FeaturePlot(seurat.integrated, reduction = "umap", features="energy", split.by="sample", ncol=min(n_samples, 4), pt.size=2) + ggtitle("- Energy")
png(file=outfile6, width=plot_width, height=750)
print(g)
dev.off()

# Plot: energy contour
df <- as.data.frame(seurat.integrated@reductions$umap@cell.embeddings)
colnames(df) <- c("UMAP_1", "UMAP_2")
df$neg_energy <- seurat.integrated$energy

g <- ggplot(df, aes(x = UMAP_1, y = UMAP_2)) +
  geom_point(color = "grey60", size = 0.5, alpha = 0.4) +
  stat_density_2d(aes(fill = ..density.., z = neg_energy),
                  geom = "tile",
                  contour = FALSE,
                  n = 100,
                  alpha = 0.45) +
  scale_fill_gradient(low = "white", high = "#FFB300") +
  geom_density_2d(aes(z = neg_energy),
                  color = "#CC8800",
                  alpha = 0.7,
                  size = 0.4) +
  theme_void() +
  theme(legend.position = 'none',
        plot.background = element_rect(fill = "white", color = NA))

ggsave(file=outfile7, g, dpi=200, width=6, height=6)
