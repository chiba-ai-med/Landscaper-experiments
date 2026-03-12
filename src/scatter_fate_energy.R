source("src/Functions.R")
library(igraph)

# Parameter
args <- commandArgs(trailingOnly=TRUE)
infile1 <- args[1]  # seurat_landscaper.RData (for major_group, E)
infile2 <- args[2]  # P_metropolis_fp.RData or P_glauber_fp.RData
infile3 <- args[3]  # igraph.RData
infile4 <- args[4]  # energy_ylim.tsv
outfile <- args[5]  # FINISH

# Loading
load(infile1)
E <- seurat.integrated@misc$landscaper$E
Group <- seurat.integrated@misc$landscaper$major_group
load(infile2)  # res_m or res_g, H_m or H_g, argmax_m or argmax_g
load(infile3)  # g (igraph object)

# Read SubGraph.tsv from same directory as igraph.RData
landscaper_dir <- dirname(infile3)
subgraph_file <- file.path(landscaper_dir, "SubGraph.tsv")
SubGraph <- unlist(read.table(subgraph_file, header = FALSE))

# Detect which dynamics (metropolis or glauber)
if (exists("res_m")) {
  res <- res_m
  argmax <- argmax_m
  dyn_name <- "metropolis"
} else if (exists("res_g")) {
  res <- res_g
  argmax <- argmax_g
  dyn_name <- "glauber"
} else {
  stop("res_m or res_g not found in RData")
}

# Y-axis range (per-PC, computed from all samples)
ylim_df <- read.table(infile4, header=TRUE, sep="\t")
YLIM <- c(ylim_df$min, ylim_df$max)

# Output directory
outdir <- dirname(outfile)
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# Basin and transient indices
Basin <- res$absorbing
transient <- res$transient
n_states <- length(E)

# Major group per state (last column of Group)
cell_type_raw <- Group[, ncol(Group)]
cell_type_raw[is.na(cell_type_raw)] <- "Uncharacterized"
cell_type <- sub("\\.[0-9]+$", "", cell_type_raw)

# Get celltype_colors from Seurat metadata
meta <- seurat.integrated@meta.data
if ("celltype_colors" %in% colnames(meta) && "celltype" %in% colnames(meta)) {
  # Build color map from metadata
  ct_color_map <- setNames(
    meta$celltype_colors[!duplicated(meta$celltype)],
    unique(meta$celltype)
  )
} else {
  # Fallback colors
  ct_color_map <- c(
    "Aboral_ectoderm" = "#E41A1C",
    "Oral_ectoderm" = "#377EB8",
    "Ciliary_band" = "#4DAF4A",
    "Neurons" = "#984EA3",
    "Uncharacterized" = "#999999"
  )
}

# Add default for Uncharacterized if not present
if (!"Uncharacterized" %in% names(ct_color_map)) {
  ct_color_map["Uncharacterized"] <- "#999999"
}

# Set factor levels
ct_levels <- unique(cell_type)

# Extract edges from igraph (observed transitions)
edge_list <- as_edgelist(g, names = FALSE)
edge_df_all <- data.frame(from = edge_list[, 1], to = edge_list[, 2])

# Get basin labels
nb <- ncol(res$F)
basin_labels <- if (!is.null(colnames(res$F))) colnames(res$F) else as.character(Basin)
basin_labels <- sub("^V", "", basin_labels)

# Get SubGraph ID for each Basin
basin_subgraph <- SubGraph[Basin]

# Create plots for each basin
for (i in seq_len(nb)) {
  target_basin <- Basin[i]
  target_subgraph <- basin_subgraph[i]

  # Filter states belonging to this basin's subgraph
  states_in_subgraph <- which(SubGraph == target_subgraph)

  # Fate probability for basin i (transient states only)
  fp_transient <- as.numeric(res$F[, i])

  # Build data frame for states in this subgraph
  df <- data.frame(
    state = seq_len(n_states),
    energy = E,
    cell_type = factor(cell_type, levels = ct_levels),
    stringsAsFactors = FALSE
  )

  # Add fate probability
  df$fate_prob <- NA
  df$fate_prob[transient] <- fp_transient
  df$fate_prob[Basin] <- ifelse(Basin == target_basin, 1.0, 0.0)

  # Mark basin vs transient
  df$is_basin <- df$state %in% Basin

  # Filter to states in this subgraph
  df_filtered <- df[df$state %in% states_in_subgraph, ]

  # Filter edges to only those within this subgraph
  edge_df <- edge_df_all[
    edge_df_all$from %in% states_in_subgraph &
    edge_df_all$to %in% states_in_subgraph,
  ]

  # Add coordinates for edges
  edge_df$x <- df_filtered$fate_prob[match(edge_df$from, df_filtered$state)]
  edge_df$y <- df_filtered$energy[match(edge_df$from, df_filtered$state)]
  edge_df$xend <- df_filtered$fate_prob[match(edge_df$to, df_filtered$state)]
  edge_df$yend <- df_filtered$energy[match(edge_df$to, df_filtered$state)]
  edge_df <- edge_df[complete.cases(edge_df), ]

  cat("Basin", basin_labels[i], "(SubGraph", target_subgraph, "): states =",
      nrow(df_filtered), ", edges =", nrow(edge_df), "\n")

  # Create plot
  p <- ggplot(df_filtered, aes(x = fate_prob, y = energy)) +
    geom_segment(
      data = edge_df,
      aes(x = x, y = y, xend = xend, yend = yend),
      color = "black",
      alpha = 0.5,
      linewidth = 1.2
    ) +
    geom_point(
      aes(color = cell_type, size = is_basin, shape = is_basin),
      alpha = 0.8
    ) +
    scale_color_manual(
      values = ct_color_map,
      name = "Cell Type",
      drop = FALSE
    ) +
    scale_size_manual(
      values = c("FALSE" = 12, "TRUE" = 24),
      guide = "none"
    ) +
    scale_shape_manual(
      values = c("FALSE" = 16, "TRUE" = 18),
      guide = "none"
    ) +
    coord_cartesian(ylim = YLIM) +
    labs(
      x = paste0("Fate Probability (Basin ", basin_labels[i], ")"),
      y = "Energy"
    ) +
    theme_classic(base_size = 28) +
    theme(
      legend.position = "right",
      plot.title = element_text(hjust = 0.5)
    )

  # Save
  fname <- file.path(outdir, paste0("scatter_fate_energy_basin_", basin_labels[i], ".png"))
  ggsave(fname, plot = p,  width = 20, height = 8, dpi = 150)
}

# Create combined plot showing ALL states (using all igraph edges)
df_all <- data.frame(
  state = seq_len(n_states),
  energy = E,
  cell_type = factor(cell_type, levels = ct_levels),
  is_basin = seq_len(n_states) %in% Basin,
  stringsAsFactors = FALSE
)

# Assign dominant basin to each state
df_all$dominant_basin <- NA
df_all$dominant_basin[transient] <- basin_labels[argmax]
df_all$dominant_basin[Basin] <- basin_labels[match(Basin, Basin)]

# Add max fate probability (how committed to dominant basin)
max_fp <- apply(res$F, 1, max)
df_all$max_fp <- NA
df_all$max_fp[transient] <- max_fp
df_all$max_fp[Basin] <- 1.0

# Build edges for all states (using original igraph edges)
edge_all <- edge_df_all
edge_all$x <- df_all$max_fp[edge_all$from]
edge_all$y <- df_all$energy[edge_all$from]
edge_all$xend <- df_all$max_fp[edge_all$to]
edge_all$yend <- df_all$energy[edge_all$to]
edge_all <- edge_all[complete.cases(edge_all), ]

# Color by dominant basin
basin_colors <- setNames(
  scales::hue_pal()(length(basin_labels)),
  basin_labels
)

p_combined <- ggplot(df_all, aes(x = max_fp, y = energy)) +
  geom_segment(
    data = edge_all,
    aes(x = x, y = y, xend = xend, yend = yend),
    color = "black",
    alpha = 0.3,
    linewidth = 0.8
  ) +
  geom_point(
    aes(color = dominant_basin, size = is_basin, shape = is_basin),
    alpha = 0.8
  ) +
  scale_color_manual(
    values = basin_colors,
    name = "Dominant Basin",
    na.value = "gray50"
  ) +
  scale_size_manual(
    values = c("FALSE" = 10, "TRUE" = 20),
    guide = "none"
  ) +
  scale_shape_manual(
    values = c("FALSE" = 16, "TRUE" = 18),
    guide = "none"
  ) +
  coord_cartesian(ylim = YLIM) +
  labs(
    x = "Max Fate Probability",
    y = "Energy"
  ) +
  theme_classic(base_size = 28) +
  theme(legend.position = "right")

fname_combined <- file.path(outdir, "scatter_fate_energy_combined.png")
ggsave(fname_combined, plot = p_combined, width = 20, height = 7, dpi = 150)

# Also create version colored by cell type with facets by dominant basin
# Filter edges: only show edges where both endpoints belong to the same dominant basin
edge_facet <- edge_all
edge_facet$dom_from <- df_all$dominant_basin[edge_facet$from]
edge_facet$dom_to <- df_all$dominant_basin[edge_facet$to]
edge_facet <- edge_facet[!is.na(edge_facet$dom_from) & !is.na(edge_facet$dom_to) &
                         edge_facet$dom_from == edge_facet$dom_to, ]
edge_facet$dominant_basin <- edge_facet$dom_from

p_celltype <- ggplot(df_all, aes(x = max_fp, y = energy)) +
  geom_segment(
    data = edge_facet,
    aes(x = x, y = y, xend = xend, yend = yend),
    color = "black",
    alpha = 0.3,
    linewidth = 1.5
  ) +
  geom_point(
    aes(color = cell_type, size = is_basin, shape = is_basin),
    alpha = 0.8
  ) +
  scale_color_manual(
    values = ct_color_map,
    name = "Cell Type",
    drop = FALSE
  ) +
  scale_size_manual(
    values = c("FALSE" = 10, "TRUE" = 20),
    guide = "none"
  ) +
  scale_shape_manual(
    values = c("FALSE" = 16, "TRUE" = 18),
    guide = "none"
  ) +
  coord_cartesian(ylim = YLIM) +
  facet_wrap(~ dominant_basin, ncol = 2) +
  labs(
    x = "Max Fate Probability",
    y = "Energy"
  ) +
  theme_classic(base_size = 24) +
  theme(
    legend.position = "right",
    strip.background = element_rect(fill = "gray95"),
    strip.text = element_text(face = "bold")
  )

fname_celltype <- file.path(outdir, "scatter_fate_energy_by_basin.png")
ggsave(fname_celltype, plot = p_celltype, width = 20, height = 7, dpi = 150)

# FINISH
file.create(outfile)
