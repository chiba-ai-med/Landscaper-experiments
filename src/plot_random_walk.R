source("src/Functions.R")

# Parameter
args <- commandArgs(trailingOnly=TRUE)
infile1  <- args[1]  # integrated seurat_landscaper.RData (for major_group labels)
infile2  <- args[2]  # random_walk list output (.RData)
infile3  <- args[3]  # current sample seurat_landscaper.RData (for Basin)

# outfiles: 10枚（5 celltype × basin/non）
outfile_ab_basin   <- args[4]
outfile_or_basin   <- args[5]
outfile_ci_basin   <- args[6]
outfile_an_basin   <- args[7]
outfile_ne_basin   <- args[8]
outfile_ab_non     <- args[9]
outfile_or_non     <- args[10]
outfile_ci_non     <- args[11]
outfile_an_non     <- args[12]
outfile_ne_non     <- args[13]

# ----------------------------
# Helpers
# ----------------------------
save_blank_png <- function(outfile, width = 8, height = 6, dpi = 300) {
  png(outfile, width = width, height = height, units = "in", res = dpi)
  par(mar = c(0, 0, 0, 0))
  plot.new()
  dev.off()
}

# ----------------------------
# Load from current sample seurat_landscaper.RData
# Use cell-level celltype to determine dominant type per state
# (major_group comes from Landscaper's hierarchical clustering and
#  may not reflect the actual cell composition of each state)
# ----------------------------
sample_env <- new.env()
load(infile3, envir = sample_env)
seurat_obj <- sample_env$seurat.integrated
basin_rows <- seurat_obj@misc$landscaper$Basin
n_states <- length(seurat_obj@misc$landscaper$E)

# Build per-state cell type membership from cell-level metadata
cell_celltype <- as.character(seurat_obj@meta.data$celltype)
cell_state_idx <- seurat_obj@meta.data$state_idx

# state_celltypes[[s]] = character vector of cell types present in state s
# (used for start state selection: a state "contains" a cell type if any cell belongs to it)
state_celltypes <- vector("list", n_states)
for (s in seq_len(n_states)) {
  ct <- cell_celltype[cell_state_idx == s]
  ct <- ct[!is.na(ct)]
  state_celltypes[[s]] <- unique(ct)
}

rm(seurat_obj, sample_env)

row_id <- seq_len(n_states)

basin_rows <- basin_rows[!is.na(basin_rows)]
basin_rows <- unique(basin_rows)
# keep only valid indices
basin_rows <- basin_rows[basin_rows >= 1 & basin_rows <= n_states]

is_basin_row <- row_id %in% basin_rows

# ----------------------------
# Load random-walk probability list
# ----------------------------
load(infile2)
if (exists("P_g_list")) {
  P_list <- P_g_list
} else if (exists("P_m_list")) {
  P_list <- P_m_list
} else {
  stop("P_g_list または P_m_list が見つかりません")
}

if (!is.list(P_list) || length(P_list) == 0) stop("P_list is empty or not a list")

# Steps (0..n)
steps <- 0:(length(P_list) - 1)

# Destination grouping
dest_levels <- c("Aboral_ectoderm", "Oral_ectoderm", "Ciliary_band", "Anus", "Neurons", "Other")
focus_types <- c("Aboral_ectoderm", "Oral_ectoderm", "Ciliary_band", "Anus", "Neurons")

# State x celltype proportion matrix
# Each row sums to 1.0 and represents the celltype composition of that state
state_ct_prop <- matrix(0, nrow = n_states, ncol = length(dest_levels),
                        dimnames = list(NULL, dest_levels))
for (s in seq_len(n_states)) {
  ct <- cell_celltype[cell_state_idx == s]
  ct <- ct[!is.na(ct)]
  if (length(ct) > 0) {
    ct2 <- ifelse(ct %in% focus_types, ct, "Other")
    tab <- table(factor(ct2, levels = dest_levels))
    state_ct_prop[s, ] <- as.numeric(tab) / sum(tab)
  }
}

# Colors (destination colors)
pal <- c(
  "Aboral_ectoderm" = "#008080",
  "Oral_ectoderm"   = "#FFFF00",
  "Ciliary_band"    = "#00FF26",
  "Anus"            = "#3F3F7F",
  "Neurons"         = "#FF00FF",
  "Other"           = "grey70"
)

plot_one_start <- function(start_label, basin_flag, outfile) {
  # Find states that contain cells of start_label type
  has_celltype <- vapply(state_celltypes, function(ct) start_label %in% ct, logical(1))
  idx_start <- which(has_celltype & is_basin_row == basin_flag)

  # If start rows not found -> blank
  if (length(idx_start) == 0) {
    save_blank_png(outfile)
    return(invisible(NULL))
  }

  # Weight starting states by the number of start_label cells they contain
  start_weights <- vapply(idx_start, function(s) {
    sum(cell_celltype[cell_state_idx == s] == start_label, na.rm = TRUE)
  }, numeric(1))
  start_weights <- start_weights / sum(start_weights)

  df_list <- vector("list", length(P_list))

  for (k in seq_along(P_list)) {
    Pk <- P_list[[k]]
    if (!is.matrix(Pk)) stop(sprintf("P_list[[%d]] is not a matrix", k))

    if (steps[k] == 0) {
      # Step 0: cells have not moved yet, so 100% is the starting cell type
      v <- setNames(rep(0, length(dest_levels)), dest_levels)
      v_label <- if (start_label %in% dest_levels) start_label else "Other"
      v[v_label] <- 1
    } else {
      # Step 1+: weighted average across starting states, then map to cell types
      P_dest_rows <- Pk[idx_start, , drop = FALSE]
      p_dest <- as.numeric(start_weights %*% P_dest_rows)

      # Distribute each state's probability proportionally across cell types
      v <- as.numeric(p_dest %*% state_ct_prop)
      names(v) <- dest_levels
    }

    df_list[[k]] <- tibble(
      step = steps[k],
      Aboral_ectoderm = v["Aboral_ectoderm"],
      Oral_ectoderm   = v["Oral_ectoderm"],
      Ciliary_band    = v["Ciliary_band"],
      Anus            = v["Anus"],
      Neurons         = v["Neurons"],
      Other           = v["Other"]
    )
  }

  df <- bind_rows(df_list) %>%
    pivot_longer(-step, names_to = "dest", values_to = "prob") %>%
    mutate(dest = factor(dest, levels = dest_levels))

  # If everything is zero -> blank
  if (all(df$prob == 0)) {
    save_blank_png(outfile)
    return(invisible(NULL))
  }

  p <- ggplot(df, aes(x = step, y = prob, color = dest)) +
    geom_line(linewidth = 2.0) +
    geom_point(size = 2.5) +
    scale_color_manual(values = pal, guide = "none") +
    scale_x_continuous(breaks = steps) +
    coord_cartesian(ylim = c(0, 1)) +
    labs(x = "Step", y = "Probability") +
    theme_classic(base_size = 12) +
    theme(text = element_text(size = 36))

  ggsave(outfile, p, width = 8, height = 6, dpi = 300)
}

# ----------------------------
# Run (5 starts × basin/nonbasin)
# ----------------------------
plot_one_start("Aboral_ectoderm", TRUE,  outfile_ab_basin)
plot_one_start("Oral_ectoderm",   TRUE,  outfile_or_basin)
plot_one_start("Ciliary_band",    TRUE,  outfile_ci_basin)
plot_one_start("Anus",            TRUE,  outfile_an_basin)
plot_one_start("Neurons",         TRUE,  outfile_ne_basin)

plot_one_start("Aboral_ectoderm", FALSE, outfile_ab_non)
plot_one_start("Oral_ectoderm",   FALSE, outfile_or_non)
plot_one_start("Ciliary_band",    FALSE, outfile_ci_non)
plot_one_start("Anus",            FALSE, outfile_an_non)
plot_one_start("Neurons",         FALSE, outfile_ne_non)