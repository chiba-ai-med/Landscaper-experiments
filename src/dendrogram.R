source("src/Functions.R")

args <- commandArgs(trailingOnly = TRUE)
infile1 <- args[1]  # seurat_landscaper.RData (for E, Basin, EnergyBarrier, Allstates_major_group)
infile2 <- args[2]  # TippingState.tsv
outfile <- args[3]

# Loading
load(infile1)
E <- seurat.integrated@misc$landscaper$E
Basin <- seurat.integrated@misc$landscaper$Basin
EnergyBarrier <- seurat.integrated@misc$landscaper$EnergyBarrier
Allstates_mg <- seurat.integrated@misc$landscaper$Allstates_major_group
TippingState <- as.matrix(read.table(infile2, header = FALSE))

n_Basin <- length(Basin)
basin_E <- E[Basin]

# Edge case: single basin
if (n_Basin <= 1) {
  hc <- NULL
  dend <- NULL
  node_order <- NULL
  df_xy <- NULL
  subtrees <- NULL
  leaves <- NULL
  paths <- NULL
  dg_skeleton <- NULL
  save(hc, dend, node_order, df_xy, subtrees, leaves, paths, dg_skeleton, file = outfile)
  quit(save = "no")
}

# Hierarchical clustering (single linkage on energy barrier)
d <- as.dist(EnergyBarrier)
hc <- hclust(d, method = "single")
dend <- as.dendrogram(hc)
node_order <- hc$order

# --- Build coordinate system ---
# Leaf x positions (in dendrogram order)
leaf_x <- numeric(n_Basin)
for (i in seq_len(n_Basin)) {
  leaf_x[node_order[i]] <- i
}

# Helper: collect all original leaf indices under a merge node
get_leaves <- function(node) {
  if (node < 0) return(-node)
  c(get_leaves(hc$merge[node, 1]), get_leaves(hc$merge[node, 2]))
}

# Internal node x positions (midpoint of children)
n_merge <- nrow(hc$merge)
internal_x <- numeric(n_merge)
internal_y <- hc$height

for (k in seq_len(n_merge)) {
  cx <- numeric(2)
  for (j in 1:2) {
    child <- hc$merge[k, j]
    if (child < 0) {
      cx[j] <- leaf_x[-child]
    } else {
      cx[j] <- internal_x[child]
    }
  }
  internal_x[k] <- mean(cx)
}

# --- Build dg_skeleton ---
# Each row: parent (x,y) -> child (xend,yend), state = basin label for leaves
skeleton_rows <- vector("list", 2 * n_merge)
idx <- 0
for (k in seq_len(n_merge)) {
  px <- internal_x[k]
  py <- internal_y[k]

  for (j in 1:2) {
    child <- hc$merge[k, j]
    if (child < 0) {
      orig <- -child
      cx <- leaf_x[orig]
      cy <- basin_E[orig]
      state_label <- as.character(Basin[orig])
    } else {
      cx <- internal_x[child]
      cy <- internal_y[child]
      state_label <- NA_character_
    }
    idx <- idx + 1
    skeleton_rows[[idx]] <- data.frame(
      x = px, y = py, xend = cx, yend = cy,
      state = state_label, stringsAsFactors = FALSE
    )
  }
}
dg_skeleton <- do.call(rbind, skeleton_rows)

# --- Build df_xy ---
df_rows <- vector("list", n_Basin + n_merge)

# Leaves (basins)
for (i in seq_len(n_Basin)) {
  df_rows[[i]] <- data.frame(
    x = leaf_x[i], y = basin_E[i],
    state = "basin",
    tipping_label = as.character(Basin[i]),
    stringsAsFactors = FALSE
  )
}

# Internal nodes (transitions)
for (k in seq_len(n_merge)) {
  # Find tipping state for this merge:
  # the pair across the two child subtrees with the minimum barrier
  leaves1 <- get_leaves(hc$merge[k, 1])
  leaves2 <- get_leaves(hc$merge[k, 2])

  tip_label <- ""
  min_barrier <- Inf
  for (l1 in leaves1) {
    for (l2 in leaves2) {
      if (!is.na(TippingState[l1, l2]) && !is.na(EnergyBarrier[l1, l2])) {
        if (EnergyBarrier[l1, l2] < min_barrier) {
          min_barrier <- EnergyBarrier[l1, l2]
          tip_label <- as.character(TippingState[l1, l2])
        }
      }
    }
  }

  df_rows[[n_Basin + k]] <- data.frame(
    x = internal_x[k], y = internal_y[k],
    state = "transition",
    tipping_label = tip_label,
    stringsAsFactors = FALSE
  )
}
df_xy <- do.call(rbind, df_rows)

# Additional objects
subtrees <- NULL
leaves <- data.frame(
  basin_idx = seq_len(n_Basin),
  basin_state = Basin,
  x = leaf_x,
  y = basin_E,
  stringsAsFactors = FALSE
)
paths <- NULL

save(hc, dend, node_order, df_xy, subtrees, leaves, paths, dg_skeleton, file = outfile)
