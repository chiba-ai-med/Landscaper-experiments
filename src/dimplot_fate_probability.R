source("src/Functions.R")

# Parameter
args <- commandArgs(trailingOnly=TRUE)
infile1 <- args[1]  # seurat_landscaper.RData (for state_idx, major_group, Seurat)
infile2 <- args[2]  # P_metropolis_fp.RData
infile3 <- args[3]  # P_glauber_fp.RData
outfile1 <- args[4]
outfile2 <- args[5]

#----------------------------------
# 出力ディレクトリ
#----------------------------------
outdir_m <- dirname(outfile1)
outdir_g <- dirname(outfile2)
dir.create(outdir_m, showWarnings = FALSE, recursive = TRUE)
dir.create(outdir_g, showWarnings = FALSE, recursive = TRUE)

#----------------------------------
# Loading
#----------------------------------
load(infile1)  # seurat.integrated
Group <- seurat.integrated@misc$landscaper$major_group
state_id <- seurat.integrated$state_idx
load(infile2)  # res_m, H_m, argmax_m, res_m$absorbing
load(infile3)  # res_g, H_g, argmax_g, res_g$absorbing

N_state <- nrow(Group)
N_cell  <- ncol(seurat.integrated)

#----------------------------------
# Basin / transient と F/H の状態レベルへの埋め戻し
#----------------------------------
Basin_m <- res_m$absorbing  # 例: c(144,176,806,607)
Basin_g <- res_g$absorbing
stopifnot(all(sort(Basin_m) == sort(Basin_g)))

transient <- sort(setdiff(seq_len(N_state), Basin_m))

stopifnot(nrow(res_m$F) == length(transient))
stopifnot(nrow(res_g$F) == length(transient))

nb_m <- ncol(res_m$F)
nb_g <- ncol(res_g$F)

# Metropolis: 状態レベルの Fate Probability (1024 x nb_m)
F_state_m <- matrix(0, nrow = N_state, ncol = nb_m)
rownames(F_state_m) <- seq_len(N_state)
F_state_m[transient, ] <- as.matrix(res_m$F)

# Basin 行：自分の Basin は 1、それ以外の Basin は 0
for (j in seq_along(Basin_m)) {
    b <- Basin_m[j]
    F_state_m[b, ] <- 0
    F_state_m[b, j] <- 1
}

# Metropolis: 状態レベル Entropy (長さ 1024)
H_state_m <- numeric(N_state)
H_state_m[transient] <- H_m
H_state_m[Basin_m]   <- 0

# Glauber: 状態レベルの Fate Probability / Entropy
F_state_g <- matrix(0, nrow = N_state, ncol = nb_g)
rownames(F_state_g) <- seq_len(N_state)
F_state_g[transient, ] <- as.matrix(res_g$F)
for (j in seq_along(Basin_g)) {
    b <- Basin_g[j]
    F_state_g[b, ] <- 0
    F_state_g[b, j] <- 1
}
H_state_g <- numeric(N_state)
H_state_g[transient] <- H_g
H_state_g[Basin_g]   <- 0

#----------------------------------
# 状態レベル → 細胞レベルへ写像
#----------------------------------
fate_cell_m <- matrix(NA_real_, nrow = N_cell, ncol = nb_m)
fate_cell_g <- matrix(NA_real_, nrow = N_cell, ncol = nb_g)
H_cell_m    <- rep(NA_real_, N_cell)
H_cell_g    <- rep(NA_real_, N_cell)

valid <- !is.na(state_id)
fate_cell_m[valid, ] <- F_state_m[state_id[valid], , drop = FALSE]
fate_cell_g[valid, ] <- F_state_g[state_id[valid], , drop = FALSE]
H_cell_m[valid]      <- H_state_m[state_id[valid]]
H_cell_g[valid]      <- H_state_g[state_id[valid]]

#----------------------------------
# Seurat meta.data に書き込み
#----------------------------------
meta <- seurat.integrated@meta.data

basin_labels_m <- as.character(Basin_m)
basin_labels_g <- as.character(Basin_g)

# Metropolis Fate Probability
for (j in seq_len(nb_m)) {
    nm <- paste0("fp_m_", basin_labels_m[j])
    meta[[nm]] <- fate_cell_m[, j]
}
meta[["entropy_m"]] <- H_cell_m

# Glauber Fate Probability
for (j in seq_len(nb_g)) {
    nm <- paste0("fp_g_", basin_labels_g[j])
    meta[[nm]] <- fate_cell_g[, j]
}
meta[["entropy_g"]] <- H_cell_g

seurat.integrated@meta.data <- meta

#----------------------------------
# プロット用の inferno カラー & テーマ
#----------------------------------
inferno_scale <- scale_color_viridis_c(option = "inferno")

clean_theme <- theme_bw() +
    theme(
        panel.grid   = element_blank(),
        panel.border = element_blank(),
        axis.text    = element_blank(),
        axis.ticks   = element_blank(),
        axis.title   = element_blank(),
        plot.title   = element_blank()
    )

#----------------------------------
# Metropolis: Basinごとの Fate Probability UMAP + Entropy UMAP
#----------------------------------
for (j in seq_len(nb_m)) {
    feat <- paste0("fp_m_", basin_labels_m[j])

    plt_list <- FeaturePlot(
        seurat.integrated,
        features  = feat,
        pt.size = 2,
        reduction = "umap",
        combine   = FALSE
    )
    p <- plt_list[[1]] + inferno_scale + clean_theme

    fname <- file.path(outdir_m, paste0("Basin_", basin_labels_m[j], "_Metropolis.png"))
    ggsave(fname, plot = p, width = 6, height = 6, dpi = 300)
}

plt_list_ent_m <- FeaturePlot(
    seurat.integrated,
    features  = "entropy_m",
    pt.size = 2,
    reduction = "umap",
    combine   = FALSE
)
p_ent_m <- plt_list_ent_m[[1]] + inferno_scale + clean_theme

ggsave(file.path(outdir_m, "Entropy_Metropolis.png"),
       plot = p_ent_m, width = 6, height = 6, dpi = 300)

#----------------------------------
# Glauber: Basinごとの Fate Probability UMAP + Entropy UMAP
#----------------------------------
for (j in seq_len(nb_g)) {
    feat <- paste0("fp_g_", basin_labels_g[j])

    plt_list <- FeaturePlot(
        seurat.integrated,
        features  = feat,
        pt.size = 2,
        reduction = "umap",
        combine   = FALSE
    )
    p <- plt_list[[1]] + inferno_scale + clean_theme

    fname <- file.path(outdir_g, paste0("Basin_", basin_labels_g[j], "_Glauber.png"))
    ggsave(fname, plot = p, width = 6, height = 6, dpi = 300)
}

plt_list_ent_g <- FeaturePlot(
    seurat.integrated,
    features  = "entropy_g",
    pt.size = 2,
    reduction = "umap",
    combine   = FALSE
)
p_ent_g <- plt_list_ent_g[[1]] + inferno_scale + clean_theme

ggsave(file.path(outdir_g, "Entropy_Glauber.png"),
       plot = p_ent_g, width = 6, height = 6, dpi = 300)

#----------------------------------
# FINISH フラグ
#----------------------------------
file.create(outfile1)
file.create(outfile2)