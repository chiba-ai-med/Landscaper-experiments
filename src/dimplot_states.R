source("src/Functions.R")

args <- commandArgs(trailingOnly = TRUE)
input_rdata <- args[1]
outfile1 <- args[2]
outfile2 <- args[3]

# Load enriched Seurat
load(input_rdata)

# Plot
g <- DimPlot(seurat.integrated, reduction = "umap", group.by="state_idx", label=FALSE, pt.size=2, label.size=6) + NoLegend() +
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
png(file=outfile1, width=600, height=600)
print(g)
dev.off()

# Plot (split by sample)
g <- DimPlot(seurat.integrated, reduction = "umap", group.by="state_idx", split.by="sample", label=FALSE, pt.size=2, label.size=6) + NoLegend() +
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
n_samples <- length(unique(seurat.integrated$sample))
plot_width <- max(500 * n_samples, 2000)
png(file=outfile2, width=plot_width, height=650)
print(g)
dev.off()
