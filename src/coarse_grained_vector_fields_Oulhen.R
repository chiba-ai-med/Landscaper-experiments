source("src/Functions.R")

args <- commandArgs(trailingOnly = TRUE)
infile1 <- args[1]
infile2 <- args[2]
infile3 <- args[3]
infile4 <- args[4]
infile5 <- args[5]
infile6 <- args[6]
outfile1 <- args[7]
outfile2 <- args[8]
# infile1 <- 'output/hpbase/integrated/seurat_annotated_landscaper.RData'
# infile2 <- 'output/hpbase/cont/seurat_annotated_landscaper.RData'
# infile3 <- 'output/hpbase/cont/sbmfcv/BIN_DATA.tsv'
# infile4 <- 'plot/hpbase/cont/Landscaper/Allstates.tsv'
# infile5 <- 'plot/hpbase/cont/P_metropolis.tsv'
# infile6 <- 'plot/hpbase/cont/P_glauber.tsv'

# Load
load(infile1)
umap <- Embeddings(seurat.integrated, "umap")[, 1:2, drop = FALSE]
load(infile2)
umap <- umap[colnames(seurat.integrated), ]

BIN <- as.matrix(read.table(infile3, header=FALSE))
Allstates <- as.matrix(read.table(infile4, header=FALSE))
P_m <- as.matrix(read.table(infile5, header=FALSE))
P_g <- as.matrix(read.table(infile6, header=FALSE))

# 2) BIN_DATA を行列化（各行が {-1,1} パターン）
#    上の固定行で unlist() されているので、必要に応じて形を直す
Sdim <- ncol(Allstates)
P <- nrow(Allstates)

# 3) 細胞→状態の対応（state_idx: 1..P）
#    文字列キーで Allstates の行と厳密一致
key <- function(v) paste0(v, collapse = "|")
keys_all   <- apply(Allstates, 1, key)
keys_cell  <- apply(BIN,        1, key)
idx_map    <- setNames(seq_len(P), keys_all)
state_idx  <- unname(idx_map[keys_cell])

# 4) 状態ごとの UMAP 重心 mu（P × 2）
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
