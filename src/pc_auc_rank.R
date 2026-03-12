#!/usr/bin/env Rscript
# pc_auc_rank.R
# Compute one-vs-rest logistic regression AUC for each PC axis and celltype
#
# Usage:
#   Rscript src/pc_auc_rank.R <pca_file> <group_file> <out_ovr> <out_max> <kfold> <seed>
#
# Input:
#   pca_file: space-separated PCA scores (no header), one row per cell
#   group_file: celltype labels (no header), one row per cell

# AUC calculation without pROC (Wilcoxon-Mann-Whitney statistic)
calc_auc <- function(y_true, y_prob) {
  pos_idx <- which(y_true == 1)
  neg_idx <- which(y_true == 0)
  n_pos <- length(pos_idx)
  n_neg <- length(neg_idx)

  if (n_pos == 0 || n_neg == 0) return(NA_real_)

  pos_probs <- y_prob[pos_idx]
  neg_probs <- y_prob[neg_idx]

  # Count pairs where pos > neg, plus 0.5 for ties
  sum_rank <- 0
  for (p in pos_probs) {
    sum_rank <- sum_rank + sum(p > neg_probs) + 0.5 * sum(p == neg_probs)
  }

  auc <- sum_rank / (n_pos * n_neg)
  return(auc)
}

# Parse arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 6) {
  stop("Usage: Rscript pc_auc_rank.R <pca_file> <group_file> <out_ovr> <out_max> <kfold> <seed>")
}

pca_file <- args[1]
group_file <- args[2]
out_ovr <- args[3]
out_max <- args[4]
K <- as.integer(args[5])
seed <- as.integer(args[6])

set.seed(seed)

# Load data
pca_dat <- read.table(pca_file, header = FALSE, stringsAsFactors = FALSE)
group_dat <- read.table(group_file, header = FALSE, stringsAsFactors = FALSE)

# Validate
if (nrow(pca_dat) != nrow(group_dat)) {
  stop("Row count mismatch between PCA file and group file")
}

# Assign column names
n_pc <- ncol(pca_dat)
colnames(pca_dat) <- paste0("PC", seq_len(n_pc))
pc_cols <- colnames(pca_dat)

celltype <- group_dat[[1]]
celltypes <- unique(celltype)
n_cells <- nrow(pca_dat)

# Assign fold IDs
fold_ids <- sample(rep(1:K, length.out = n_cells))

# Function to compute CV AUC for one PC + one celltype (OvR)
compute_cv_auc <- function(z, y_binary, fold_ids, K) {
  preds <- rep(NA_real_, length(y_binary))

  for (k in seq_len(K)) {
    test_idx <- which(fold_ids == k)
    train_idx <- which(fold_ids != k)

    y_train <- y_binary[train_idx]

    # Skip fold if only one class in training
    if (length(unique(y_train)) < 2) next

    z_train <- z[train_idx]
    z_test <- z[test_idx]

    # Fit logistic regression
    fit <- tryCatch(
      glm(y_train ~ z_train, family = binomial()),
      warning = function(w) NULL,
      error = function(e) NULL
    )

    if (is.null(fit)) next

    # Predict
    preds[test_idx] <- tryCatch(
      predict(fit, newdata = data.frame(z_train = z_test), type = "response"),
      error = function(e) NA_real_
    )
  }

  # Compute AUC from out-of-fold predictions
  valid_idx <- which(!is.na(preds))
  y_valid <- y_binary[valid_idx]
  p_valid <- preds[valid_idx]

  # Need both classes present
  if (length(unique(y_valid)) < 2) {
    return(NA_real_)
  }

  auc_val <- tryCatch(
    calc_auc(y_valid, p_valid),
    error = function(e) NA_real_
  )

  return(auc_val)
}

# Main loop: collect results
results_ovr <- vector("list", length(pc_cols) * length(celltypes))
idx <- 1

for (pc in pc_cols) {
  z <- pca_dat[[pc]]

  for (ct in celltypes) {
    y_binary <- as.integer(celltype == ct)
    n_pos <- sum(y_binary)
    n_neg <- sum(1 - y_binary)

    auc_val <- compute_cv_auc(z, y_binary, fold_ids, K)

    results_ovr[[idx]] <- data.frame(
      pc = pc,
      celltype_pos = ct,
      auc = auc_val,
      n_pos = n_pos,
      n_neg = n_neg,
      stringsAsFactors = FALSE
    )
    idx <- idx + 1
  }
}

df_ovr <- do.call(rbind, results_ovr)

# Compute max AUC per PC
df_max <- do.call(rbind, lapply(pc_cols, function(pc) {
  sub <- df_ovr[df_ovr$pc == pc, ]
  if (all(is.na(sub$auc))) {
    data.frame(pc = pc, auc_max = NA_real_, best_celltype_pos = NA_character_,
               stringsAsFactors = FALSE)
  } else {
    best_row <- sub[which.max(sub$auc), ]
    data.frame(pc = pc, auc_max = best_row$auc, best_celltype_pos = best_row$celltype_pos,
               stringsAsFactors = FALSE)
  }
}))

# Create output directories if needed
dir.create(dirname(out_ovr), showWarnings = FALSE, recursive = TRUE)
dir.create(dirname(out_max), showWarnings = FALSE, recursive = TRUE)

# Write outputs
write.table(df_ovr, out_ovr, sep = "\t", row.names = FALSE, quote = FALSE)
write.table(df_max, out_max, sep = "\t", row.names = FALSE, quote = FALSE)

cat("Done. Wrote", nrow(df_ovr), "rows to", out_ovr, "\n")
cat("Done. Wrote", nrow(df_max), "rows to", out_max, "\n")
