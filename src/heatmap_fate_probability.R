source("src/Functions.R")

# Parameter
args <- commandArgs(trailingOnly=TRUE)
infile1 <- args[1]  # seurat_landscaper.RData (for major_group)
infile2 <- args[2]  # P_metropolis_fp.RData
infile3 <- args[3]  # P_glauber_fp.RData
outfile1 <- args[4]
outfile2 <- args[5]

# Loading
load(infile1)
Group <- seurat.integrated@misc$landscaper$major_group
load(infile2) # res_m, H_m, argmax_m
load(infile3) # res_g, H_g, argmax_g

# 出力ディレクトリ
outdir_m <- dirname(outfile1)
outdir_g <- dirname(outfile2)
dir.create(outdir_m, showWarnings = FALSE, recursive = TRUE)
dir.create(outdir_g, showWarnings = FALSE, recursive = TRUE)

# ---- Basin / transient の決定 ----
Basin <- res_m$absorbing        # 4つの Basin インデックス
transient <- sort(setdiff(seq_len(nrow(Group)), Basin))

# ---- 細胞型 one-hot（まず全1024行で作る）----
cell_type_raw <- Group[, ncol(Group)]
cell_type_raw[is.na(cell_type_raw)] <- "Unknown"
cell_type_full <- sub("\\.[0-9]+$", "", cell_type_raw)

# ★ 行順を固定
desired_order <- c(
  "Aboral_ectoderm",
  "Oral_ectoderm",
  "Ciliary_band",
  "Neurons",
  "Unknown"
)

# データに存在するものだけ
ct_levels <- intersect(desired_order, unique(cell_type_full))

X_full <- sapply(ct_levels, function(ct) as.integer(cell_type_full == ct))
X_full <- as.matrix(X_full)
colnames(X_full) <- ct_levels

# transient 行だけ抜き出す
X_type <- X_full[transient, , drop = FALSE]

# 行名
rn <- as.character(transient)
rownames(X_type)  <- rn
rownames(res_m$F) <- rn
rownames(res_g$F) <- rn

# ---- 色設定（0=暗い, 1=黄色）----
pal_binary    <- c("#777777", "#00ff00")
binary_breaks <- c(-0.5, 0.5, 1.5)

# ---- Basin 情報 ----
nb_m <- ncol(res_m$F)
nb_g <- ncol(res_g$F)
basin_labels_m <- if (!is.null(colnames(res_m$F))) colnames(res_m$F) else as.character(res_m$absorbing)
basin_labels_g <- if (!is.null(colnames(res_g$F))) colnames(res_g$F) else as.character(res_g$absorbing)
basin_labels_m <- sub("^V", "", basin_labels_m)
basin_labels_g <- sub("^V", "", basin_labels_g)

# ---- Plot (Metropolis) ----
lapply(seq_len(nb_m), function(i) {
    f_i <- as.numeric(res_m$F[, i])          # fate prob for basin i
    ord <- order(f_i, decreasing = FALSE)     # 小さい順に並べる

    # 行=細胞型, 列=細胞 にするため転置
    X_ord <- t(X_type[ord, , drop = FALSE])

    fname <- file.path(outdir_m, paste0("Basin_", basin_labels_m[i], ".png"))

    png(fname, width = 1600, height = 600)
    pheatmap(
        X_ord,
        cluster_rows   = FALSE,
        cluster_cols   = FALSE,
        show_rownames  = TRUE,    # 左に細胞型ラベル
        show_colnames  = FALSE,   # 上にセル名は出さない
        color          = pal_binary,
        breaks         = binary_breaks,
        border_color   = NA,
        legend         = FALSE
        # main は指定しない
    )
    dev.off()
})

# ---- Plot (Glauber) ----
lapply(seq_len(nb_g), function(i) {
    f_i <- as.numeric(res_g$F[, i])
    ord <- order(f_i, decreasing = TRUE)

    X_ord <- t(X_type[ord, , drop = FALSE])

    fname <- file.path(outdir_g, paste0("Basin_", basin_labels_g[i], ".png"))

    png(fname, width = 1600, height = 600)
    pheatmap(
        X_ord,
        cluster_rows   = FALSE,
        cluster_cols   = FALSE,
        show_rownames  = TRUE,
        show_colnames  = FALSE,
        color          = pal_binary,
        breaks         = binary_breaks,
        border_color   = NA,
        legend         = FALSE
    )
    dev.off()
})

# FINISH
file.create(outfile1)
file.create(outfile2)