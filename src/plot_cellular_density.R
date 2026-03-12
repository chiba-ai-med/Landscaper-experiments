source("src/Functions.R")

# Parameter
args <- commandArgs(trailingOnly = TRUE)
infile <- args[1]
outfile <- args[2]

# Loading
load(infile)

# Pre-processing
df <- as.data.frame(seurat.integrated@reductions$umap@cell.embeddings)
colnames(df) <- c("UMAP_1", "UMAP_2")

# Plot
g <- ggplot(df, aes(x = UMAP_1, y = UMAP_2)) +
  geom_point(color = "grey60", size = 0.5, alpha = 0.4) +
  stat_density_2d(aes(fill = ..density..), 
                  geom = "tile", 
                  contour = FALSE,
                  n = 100,
                  alpha = 0.6) +  # 少し透明度を上げた
  scale_fill_gradient(low = "white", high = "#0066f4ff") +  # より濃い青
  geom_density_2d(color = "#003366", alpha = 0.6, size = 0.4) +  # 濃い青の等高線
  theme_void() +
  theme(legend.position = 'none',
        plot.background = element_rect(fill = "white", color = NA))

ggsave(file=outfile, g, dpi=200, width=6, height=6)