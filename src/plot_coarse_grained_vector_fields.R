source("src/Functions.R")

args <- commandArgs(trailingOnly = TRUE)
infile1 <- args[1]  # integrated seurat_landscaper.RData (for UMAP + meta)
infile2 <- args[2]  # sample seurat_landscaper.RData (for state_idx, E)
infile3 <- args[3]  # P_metropolis_cg.RData
infile4 <- args[4]  # P_glauber_cg.RData
outfile1 <- args[5]
outfile2 <- args[6]
outfile3 <- args[7]
outfile4 <- args[8]
outfile5 <- args[9]
outfile6 <- args[10]
outfile7 <- args[11]
outfile8 <- args[12]
outfile9 <- args[13]
outfile10 <- args[14]
outfile11 <- args[15]
outfile12 <- args[16]
outfile13 <- args[17]
outfile14 <- args[18]

# Load integrated (for meta + UMAP)
int_env <- new.env()
load(infile1, envir = int_env)
meta <- int_env$seurat.integrated@meta.data
umap <- Embeddings(int_env$seurat.integrated, "umap")[, 1:2, drop = FALSE]

# Load sample (for state_idx, E)
load(infile2)
meta <- meta[colnames(seurat.integrated), ]
umap <- umap[colnames(seurat.integrated), ]
state_idx <- seurat.integrated$state_idx
E <- seurat.integrated@misc$landscaper$E
rm(int_env)

load(infile3)  # mu, have, dr_m
load(infile4)  # mu, have, dr_g

df_m <- make_arrow_df(dr_m, mu, umap, arrow_color = "#000000ff")
df_g <- make_arrow_df(dr_g, mu, umap, arrow_color = "#000000ff")

################################################
################## metropolis ##################
################################################
# none
p_bg <- base_plot(umap, mode = "none", meta = meta, E = E, germlayer_cols = NULL)

p <- p_bg +
  { if (nrow(df_m) > 0) geom_segment(data = df_m,
      aes(x = x, y = y, xend = xend, yend = yend),
      arrow = arrow(length = unit(0.35, "cm"), angle = 25, type = "open"),
      color = df_m$col, alpha = scales::rescale(df_m$mass, to = c(0.3, 1)),
      linewidth = 1.2) }

ggsave(outfile1, p, width = 7, height = 6, dpi = 300)

# germlayer
p_bg <- base_plot(umap, mode = "germlayer", meta = meta, E = E, germlayer_cols = NULL)

p <- p_bg +
  { if (nrow(df_m) > 0) geom_segment(data = df_m,
      aes(x = x, y = y, xend = xend, yend = yend),
      arrow = arrow(length = unit(0.35, "cm"), angle = 25, type = "open"),
      color = df_m$col, alpha = scales::rescale(df_m$mass, to = c(0.3, 1)),
      linewidth = 1.2) }

ggsave(outfile2, p, width = 7, height = 6, dpi = 300)

# cluster
p_bg <- base_plot(umap, mode = "cluster", meta = meta, E = E, germlayer_cols = NULL)

p <- p_bg +
  { if (nrow(df_m) > 0) geom_segment(data = df_m,
      aes(x = x, y = y, xend = xend, yend = yend),
      arrow = arrow(length = unit(0.35, "cm"), angle = 25, type = "open"),
      color = df_m$col, alpha = scales::rescale(df_m$mass, to = c(0.3, 1)),
      linewidth = 1.2) }

ggsave(outfile3, p, width = 7, height = 6, dpi = 300)

# sample
p_bg <- base_plot(umap, mode = "sample", meta = meta, E = E, germlayer_cols = NULL)

p <- p_bg +
  { if (nrow(df_m) > 0) geom_segment(data = df_m,
      aes(x = x, y = y, xend = xend, yend = yend),
      arrow = arrow(length = unit(0.35, "cm"), angle = 25, type = "open"),
      color = df_m$col, alpha = scales::rescale(df_m$mass, to = c(0.3, 1)),
      linewidth = 1.2) }

ggsave(outfile4, p, width = 7, height = 6, dpi = 300)

# celltype
p_bg <- base_plot(umap, mode = "celltype", meta = meta, E = E, germlayer_cols = NULL)

p <- p_bg +
  { if (nrow(df_m) > 0) geom_segment(data = df_m,
      aes(x = x, y = y, xend = xend, yend = yend),
      arrow = arrow(length = unit(0.35, "cm"), angle = 25, type = "open"),
      color = df_m$col, alpha = scales::rescale(df_m$mass, to = c(0.3, 1)),
      linewidth = 1.2) }

ggsave(outfile5, p, width = 7, height = 6, dpi = 300)

# state
p_bg <- base_plot(umap, mode = "state", meta = meta, E = E, germlayer_cols = NULL, state_idx=state_idx)

p <- p_bg +
  { if (nrow(df_m) > 0) geom_segment(data = df_m,
      aes(x = x, y = y, xend = xend, yend = yend),
      arrow = arrow(length = unit(0.35, "cm"), angle = 25, type = "open"),
      color = df_m$col, alpha = scales::rescale(df_m$mass, to = c(0.3, 1)),
      linewidth = 1.2) }

ggsave(outfile6, p, width = 7, height = 6, dpi = 300)

# energy
p_bg <- base_plot(umap, mode = "energy", meta = meta, E = E, germlayer_cols = NULL, state_idx=state_idx)

p <- p_bg +
  { if (nrow(df_m) > 0) geom_segment(data = df_m,
      aes(x = x, y = y, xend = xend, yend = yend),
      arrow = arrow(length = unit(0.35, "cm"), angle = 25, type = "open"),
      color = df_m$col, alpha = scales::rescale(df_m$mass, to = c(0.3, 1)),
      linewidth = 1.2) }

ggsave(outfile7, p, width = 7, height = 6, dpi = 300)


################################################
################## glauber ##################
################################################
# none
p_bg <- base_plot(umap, mode = "none", meta = meta, E = E, germlayer_cols = NULL)

p <- p_bg +
  { if (nrow(df_g) > 0) geom_segment(data = df_g,
      aes(x = x, y = y, xend = xend, yend = yend),
      arrow = arrow(length = unit(0.35, "cm"), angle = 25, type = "open"),
      color = df_g$col, alpha = scales::rescale(df_g$mass, to = c(0.3, 1)),
      linewidth = 1.2) }

ggsave(outfile8, p, width = 7, height = 6, dpi = 300)

# germlayer
p_bg <- base_plot(umap, mode = "germlayer", meta = meta, E = E, germlayer_cols = NULL)

p <- p_bg +
  { if (nrow(df_g) > 0) geom_segment(data = df_g,
      aes(x = x, y = y, xend = xend, yend = yend),
      arrow = arrow(length = unit(0.35, "cm"), angle = 25, type = "open"),
      color = df_g$col, alpha = scales::rescale(df_g$mass, to = c(0.3, 1)),
      linewidth = 1.2) }

ggsave(outfile9, p, width = 7, height = 6, dpi = 300)

# cluster
p_bg <- base_plot(umap, mode = "cluster", meta = meta, E = E, germlayer_cols = NULL)

p <- p_bg +
  { if (nrow(df_g) > 0) geom_segment(data = df_g,
      aes(x = x, y = y, xend = xend, yend = yend),
      arrow = arrow(length = unit(0.35, "cm"), angle = 25, type = "open"),
      color = df_g$col, alpha = scales::rescale(df_g$mass, to = c(0.3, 1)),
      linewidth = 1.2) }

ggsave(outfile10, p, width = 7, height = 6, dpi = 300)

# sample
p_bg <- base_plot(umap, mode = "sample", meta = meta, E = E, germlayer_cols = NULL)

p <- p_bg +
  { if (nrow(df_g) > 0) geom_segment(data = df_g,
      aes(x = x, y = y, xend = xend, yend = yend),
      arrow = arrow(length = unit(0.35, "cm"), angle = 25, type = "open"),
      color = df_g$col, alpha = scales::rescale(df_g$mass, to = c(0.3, 1)),
      linewidth = 1.2) }

ggsave(outfile11, p, width = 7, height = 6, dpi = 300)

# celltype
p_bg <- base_plot(umap, mode = "celltype", meta = meta, E = E, germlayer_cols = NULL)

p <- p_bg +
  { if (nrow(df_g) > 0) geom_segment(data = df_g,
      aes(x = x, y = y, xend = xend, yend = yend),
      arrow = arrow(length = unit(0.35, "cm"), angle = 25, type = "open"),
      color = df_g$col, alpha = scales::rescale(df_g$mass, to = c(0.3, 1)),
      linewidth = 1.2) }

ggsave(outfile12, p, width = 7, height = 6, dpi = 300)

# state
p_bg <- base_plot(umap, mode = "state", meta = meta, E = E, germlayer_cols = NULL, state_idx=state_idx)

p <- p_bg +
  { if (nrow(df_g) > 0) geom_segment(data = df_g,
      aes(x = x, y = y, xend = xend, yend = yend),
      arrow = arrow(length = unit(0.35, "cm"), angle = 25, type = "open"),
      color = df_g$col, alpha = scales::rescale(df_g$mass, to = c(0.3, 1)),
      linewidth = 1.2) }

ggsave(outfile13, p, width = 7, height = 6, dpi = 300)

# energy
p_bg <- base_plot(umap, mode = "energy", meta = meta, E = E, germlayer_cols = NULL, state_idx=state_idx)

p <- p_bg +
  { if (nrow(df_g) > 0) geom_segment(data = df_g,
      aes(x = x, y = y, xend = xend, yend = yend),
      arrow = arrow(length = unit(0.35, "cm"), angle = 25, type = "open"),
      color = df_g$col, alpha = scales::rescale(df_g$mass, to = c(0.3, 1)),
      linewidth = 1.2) }

ggsave(outfile14, p, width = 7, height = 6, dpi = 300)
