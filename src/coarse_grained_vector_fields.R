source("src/Functions.R")

args <- commandArgs(trailingOnly = TRUE)
infile1 <- args[1]  # integrated seurat_landscaper.RData (for UMAP)
infile2 <- args[2]  # sample seurat_landscaper.RData (for state_idx, Allstates)
infile3 <- args[3]  # P_metropolis.tsv
infile4 <- args[4]  # P_glauber.tsv
outfile1 <- args[5]
outfile2 <- args[6]

# Load integrated (for UMAP)
int_env <- new.env()
load(infile1, envir = int_env)
umap <- Embeddings(int_env$seurat.integrated, "umap")[, 1:2, drop = FALSE]

# Load sample (for state_idx, Allstates)
load(infile2)
umap <- umap[colnames(seurat.integrated), ]
state_idx <- seurat.integrated$state_idx
Allstates <- seurat.integrated@misc$landscaper$Allstates
rm(int_env)

P_m <- as.matrix(read.table(infile3, header=FALSE))
P_g <- as.matrix(read.table(infile4, header=FALSE))

P <- nrow(Allstates)

# 状態ごとの UMAP 重心 mu（P × 2）
# rowsum で高速集計 → 件数で割る
# counts: 各状態のセル数（長さ P）
# have: 出現している状態の論理ベクトル（長さ P）
counts <- tabulate(state_idx, nbins = P)   # 1..P の整数を仮定
have   <- counts > 0

# rowsum は「出現した状態だけ」を返す（行名=状態番号）
rs <- rowsum(umap, group = state_idx)      # デフォルト reorder=TRUE で昇順
idx_present <- as.integer(rownames(rs))    # 実際に現れた状態の添字

# P×2 に拡張した合計行列を用意して埋める
sum_xy <- matrix(0, nrow = P, ncol = 2)
sum_xy[idx_present, ] <- as.matrix(rs)

# 重心 mu（P×2）。存在しない状態は NA のままにする
mu <- matrix(NA_real_, nrow = P, ncol = 2)
mu[have, ] <- sum_xy[have, ] / counts[have]
colnames(mu) <- c("x", "y")

# 5) 要約ドリフト（行列積で一括）
dr_m <- drift_from_P(P_m, mu, counts)
dr_g <- drift_from_P(P_g, mu, counts)

# Save
save(mu, have, dr_m, file=outfile1)
save(mu, have, dr_g, file=outfile2)
