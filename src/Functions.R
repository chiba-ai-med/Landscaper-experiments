library("GSEABase")
library("GO.db")
library("Seurat")
library("sctransform")
library("readxl")
library("writexl")
library("scDblFinder")
library("scran")
library("scater")
library("SeuratWrappers")
library("monocle3")
library("velociraptor")
library("scTGIF")
library("RColorBrewer")
library("reticulate")
library("abind")
library("purrr")
library("tidyr")
library("stringr")
library("dplyr")
library("tibble")
library("magrittr")
library("ggplot2")
library("ggforce")
library("gridExtra")
library("viridis")
library("reshape2")
library("SingleCellExperiment")
library("schex")
library("hexbin")
library("metR")
library("grid")
library("akima")
library("qvalue")
library("patchwork")
library("Matrix")
library("igraph")
library("dendextend")
library("Rcpp")
library("pheatmap")

# 1) 行列P（行 stochastic）と吸収状態のインデックスから Q, R を切り出し
build_QR <- function(P, absorbing) {
  stopifnot(nrow(P) == ncol(P))
  n <- nrow(P); absorbing <- sort(unique(absorbing))
  transient <- setdiff(seq_len(n), absorbing)
  if (length(transient) == 0) stop("Transient states are empty.")
  Q <- P[transient, transient, drop = FALSE]
  R <- P[transient, absorbing, drop = FALSE]
  list(Q = as(Q, "dgCMatrix"), R = as(R, "dgCMatrix"),
       transient = transient, absorbing = absorbing)
}

# 2) 吸収確率 F を計算（F_{ij}: transient i → absorbing j）
#    逆行列は使わず (I - Q) X = R を解く
absorption_probabilities <- function(P, absorbing, reg = 0) {
  # 事前に行正規化（数値誤差対策）
  P <- as(P, "dgCMatrix")
  rs <- rowSums(P)
  if (any(rs == 0)) stop("Some rows of P sum to 0.")
  P <- Diagonal(x = 1/rs) %*% P

  qr <- build_QR(P, absorbing)
  Q <- qr$Q; R <- qr$R
  m <- nrow(Q)
  IminusQ <- Diagonal(m) - Q
  if (reg > 0) IminusQ <- IminusQ + Diagonal(m, x = reg)  # 退避的レギュラリゼーション

  # (I - Q) * F = R を各列について解く
  F <- solve(IminusQ, R, sparse = TRUE)
  # F は transient×absorbing。各行は吸収確率分布（和=1になるのが理想）
  list(F = F, transient = qr$transient, absorbing = qr$absorbing)
}

# 3) （任意）可塑性: エントロピー H_i = -∑_j F_ij log F_ij
fate_entropy <- function(F, base = exp(1)) {
  F <- as.matrix(F)
  F[F <= 0] <- 0
  rowSums(ifelse(F > 0, -F * (log(F) / log(base)), 0))
}

# 4) （任意）最大確率の運命割当（argmax）
fate_argmax <- function(F) {
  max.col(as.matrix(F), ties.method = "first")
}

make_vector_grid <- function(df, umap, nx = 80, ny = 80) {
  if (nrow(df) == 0) return(data.frame())
  dxp <- df$xend - df$x
  dyp <- df$yend - df$y
  xr <- range(umap[,1], na.rm = TRUE); yr <- range(umap[,2], na.rm = TRUE)
  xo <- seq(xr[1], xr[2], length.out = nx)
  yo <- seq(yr[1], yr[2], length.out = ny)

  ix <- akima::interp(x = df$x, y = df$y, z = dxp,
                      xo = xo, yo = yo, duplicate = "mean", linear = TRUE)
  iy <- akima::interp(x = df$x, y = df$y, z = dyp,
                      xo = xo, yo = yo, duplicate = "mean", linear = TRUE)

  grd <- expand.grid(x = ix$x, y = ix$y)
  grd$dx <- as.vector(ix$z)
  grd$dy <- as.vector(iy$z)

  # ← ここで行を落とさず、NA を 0 に置換して「規則格子」を保つ
  grd$dx[!is.finite(grd$dx)] <- 0
  grd$dy[!is.finite(grd$dy)] <- 0

  # 正規化（任意）
  sp <- sqrt(grd$dx^2 + grd$dy^2)
  smax <- max(sp, na.rm = TRUE); if (!is.finite(smax) || smax == 0) smax <- 1
  grd$dx <- grd$dx / smax
  grd$dy <- grd$dy / smax
  grd$speed <- sp / smax

  grd
}

# ---- ベースの UMAP 散布 ----
base_plot <- function(umap, mode = "none", meta = NULL, E = NULL,
                      germlayer_cols = NULL, state_idx = NULL) {
  stopifnot(is.matrix(umap), ncol(umap) >= 2)
  df <- data.frame(UMAP_1 = umap[,1], UMAP_2 = umap[,2], stringsAsFactors = FALSE)
  p  <- ggplot()

  default_disc <- function(n) scales::hue_pal()(n)
  if (is.null(germlayer_cols)) {
    germlayer_cols <- c(Ectoderm="#1b9e77", Endoderm="#d95f02", Mesoderm="#7570b3", NA_="#00000080")
  }

  if (mode == "germlayer") {
    if (is.null(meta) || !("germlayer" %in% colnames(meta)))
      stop("meta に 'germlayer' 列がありません。")
    df$grp <- factor(meta$germlayer, levels = c("Ectoderm","Endoderm","Mesoderm","NA"))
    pal <- germlayer_cols; names(pal)[names(pal)=="NA_"] <- "NA"
    p <- p + geom_point(data = df, aes(UMAP_1, UMAP_2, color = grp), size = 0.4, alpha = 0.7) +
      scale_color_manual(values = pal, name = NULL)

  } else if (mode == "cluster") {
    stopifnot(!is.null(meta))
    grp <- if ("seurat_clusters" %in% colnames(meta)) meta$seurat_clusters else Seurat::Idents(seurat.integrated)
    df$grp <- factor(as.character(grp))
    p <- p + geom_point(data = df, aes(UMAP_1, UMAP_2, color = grp), size = 0.4, alpha = 0.7) +
      scale_color_manual(values = default_disc(nlevels(df$grp)), name = "cluster")

  } else if (mode == "sample") {
    stopifnot(!is.null(meta))
    if (all(c("condition","time") %in% colnames(meta))) {
      df$grp <- factor(paste0(meta$condition, "_", meta$time))
    } else if ("sample" %in% colnames(meta)) {
      df$grp <- factor(meta$sample)
    } else {
      stop("meta に 'sample' または ('condition','time') がありません。")
    }
    p <- p + geom_point(data = df, aes(UMAP_1, UMAP_2, color = grp), size = 0.4, alpha = 0.7) +
      scale_color_manual(values = default_disc(nlevels(df$grp)), name = "sample")

} else if (mode == "celltype") {
    if (is.null(meta) || !("celltype" %in% colnames(meta)))
      stop("meta に 'celltype' 列がありません。")
    
    df$grp <- factor(meta$celltype)
    
    # celltype_colorsが存在する場合、それを使用
    if ("celltype_colors" %in% colnames(meta)) {
      # ユニークなcelltypeと対応する色を取得
      celltype_color_map <- setNames(
        meta$celltype_colors[!duplicated(meta$celltype)],
        unique(meta$celltype)
      )
      
      p <- p + geom_point(data = df, aes(UMAP_1, UMAP_2, color = grp), size = 0.4, alpha = 0.7) +
        scale_color_manual(values = celltype_color_map, name = "celltype")
    } else {
      # celltype_colorsがない場合はデフォルト
      p <- p + geom_point(data = df, aes(UMAP_1, UMAP_2, color = grp), size = 0.4, alpha = 0.7) +
        scale_color_manual(values = default_disc(nlevels(df$grp)), name = "celltype")
    }

  } else if (mode == "state") {
    if (is.null(state_idx))
      stop("state モードには state_idx が必要です（BIN×Allstatesで作ったもの）。")
    df$grp <- factor(state_idx)
    p <- p + geom_point(data = df, aes(UMAP_1, UMAP_2, color = grp), size = 0.4, alpha = 0.7) +
      scale_color_manual(values = default_disc(nlevels(df$grp)), name = "state")
  } else if (mode == "energy") {
    if (is.null(state_idx)) stop("energy モードには state_idx が必要です。")
    if (is.null(E)) stop("energy モードには E が必要です。")
    if (max(state_idx, na.rm = TRUE) > length(E))
      stop("E の長さと state_idx が不整合です。")
    df$val <- as.numeric(E[state_idx])
    p <- p + geom_point(data = df, aes(UMAP_1, UMAP_2, color = val), size = 0.4, alpha = 0.7) +
      scale_color_viridis_c(name = "Energy")

  } else {
    p <- p + geom_point(data = df, aes(UMAP_1, UMAP_2), color = "grey85", size = 0.4, alpha = 0.6)
  }
  p <- p + guides(color = "none", fill = "none") + theme(legend.position = "none")
  p + coord_equal() + theme_classic(base_size = 12) +
    theme(legend.position = "right")
}


# ---- 矢印データフレームを作る（観測状態のみ） ----
make_arrow_df <- function(dr, mu, umap, arrow_color = "steelblue") {
  P <- nrow(mu)
  if (nrow(dr$V) != P) stop("dr$V と mu の行数が一致しません")
  have <- as.logical(dr$have); have[is.na(have)] <- FALSE
  ok <- have & !is.na(mu[,1]) & !is.na(mu[,2])
  if (!any(ok)) return(data.frame())

  xr <- range(umap[,1], na.rm = TRUE); yr <- range(umap[,2], na.rm = TRUE)
  span <- max(diff(xr), diff(yr)); scale_len <- 0.10 * span

  V <- dr$V
  len <- sqrt(rowSums(V^2))
  max_len <- max(len[ok], na.rm = TRUE); if (max_len == 0) max_len <- 1
  vx <- (V[,1] / max_len) * scale_len; vy <- (V[,2] / max_len) * scale_len

  data.frame(
    x    = mu[ok, 1],
    y    = mu[ok, 2],
    xend = mu[ok, 1] + vx[ok],
    yend = mu[ok, 2] + vy[ok],
    mass = dr$mass[ok],
    col  = arrow_color
)
}

# 5) 要約ドリフト（行列積で一括）
#    Toff = P（対角0）、質量 rs = rowSums(Toff)、
#    V = Toff %*% mu - mu * cbind(rs,rs)
drift_from_P <- function(P, mu, counts = NULL) {
  Toff <- P
  diag(Toff) <- 0

  # 存在する状態フラグ（countsがあれば使う。無ければmuのNAで判定）
  if (is.null(counts)) {
    have <- rowSums(is.na(mu)) == 0
  } else {
    have <- counts > 0
  }

  # 1) 存在しない状態の座標を 0 に
  mu0 <- mu
  mu0[!have, ] <- 0

  # 2) 存在しない状態との遷移をカット（行・列をゼロ）
  if (any(!have)) {
    Toff[, !have] <- 0  # 行き先が未定義 → 合力に使わない
    Toff[!have, ] <- 0  # 出発も未定義 → ドリフトを計算しない
  }

  # 3) ドリフトと強さ
  rs <- rowSums(Toff)                 # 非対角質量（強さ）
  V  <- Toff %*% mu0 - mu0 * cbind(rs, rs)  # P×2

  list(V = V, mass = rs, Toff = Toff, have = have)
}

#################################################
# Disconnectivity graph
#################################################
# a helper function to find route for each leaf, and apply it to all leaves.
pathRoutes <- function(leaf, subtrees) {
  which(sapply(subtrees, function(x) leaf %in% x))
}

hamming1_adjacency <- function(S) {
  P <- nrow(S)
  # 各状態ベクトルを {0,1} に変換して XOR の総和で距離を取る
  # (−1,1) → (0,1)
  S01 <- (S + 1L) / 2L
  # 全ての組み合わせのハミング距離
  # 行列積を使うと速い：距離 = 異なるビット数 = 列数 − 一致数
  matches <- S01 %*% t(S01) + (1 - S01) %*% t(1 - S01)  # 一致ビット数
  dH <- ncol(S) - matches
  # ハミング距離が1の要素を1、それ以外を0
  A <- (dH == 1) * 1
  diag(A) <- 0
  A
}

sparsify_by_threshold <- function(Tmat, tau = 1e-3, mode = c("abs","quantile")) {
  stopifnot(inherits(Tmat, "dgCMatrix"))
  mode <- match.arg(mode)
  A <- Tmat

  n <- nrow(A)
  # 対角をいったん 0（あとで行和=1に戻す）
  diag(A) <- 0

  # しきい値の決定
  thr <- if (mode == "abs") {
    tau
  } else {
    xnz <- A@x
    if (!length(xnz)) 0 else {
      qs <- quantile(xnz, probs = tau, names = FALSE)
      if (is.na(qs)) 0 else qs
    }
  }

  # しきい値未満を0に
  if (thr > 0) {
    A@x[A@x < thr] <- 0
    A <- drop0(A)
  }

  # 行ごとに非対角の合計質量
  row_mass <- rowSums(A)
  # 対角に「残りの質量」を入れて行和=1に
  diag(A) <- pmax(0, 1 - row_mass)

  # 行が全ゼロの保険：自己遷移=1
  zero_rows <- which(row_mass == 0)
  if (length(zero_rows)) diag(A)[zero_rows] <- 1

  A
}

##------------------------------------------------------------------------------
## 転移行列の生成（Metropolis / Glauber）＋ 時系列制約オプション付き
##------------------------------------------------------------------------------
transition_matrix_from_E <- function(
  E,
  S,
  beta   = 1,
  kernel = c("metropolis", "glauber"),
  time_state = NULL,   # 長さPの「状態ごとの時間」(例: 36,48,72,96 の平均など)
  tol = 0              # time_state[j] < time_state[i] - tol を禁止（単位は time_state と同じ）
) {
  kernel <- match.arg(kernel)
  P     <- nrow(S)
  Sdim  <- ncol(S)

  # 状態 → 行番号の辞書（ハッシュ）
  key  <- function(v) paste0(v, collapse = "")
  keys <- apply(S, 1, key)
  idx_map <- setNames(seq_len(P), keys)

  ii <- integer(0)
  jj <- integer(0)
  xx <- numeric(0)   # triplet 形式で sparse に詰める

  for (p in seq_len(P)) {
    s <- S[p, ]

    # 1ビット反転した近傍 Sdim 個を作る
    nei <- matrix(rep(s, each = Sdim), nrow = Sdim)
    nei[cbind(seq_len(Sdim), seq_len(Sdim))] <- -nei[cbind(seq_len(Sdim), seq_len(Sdim))]

    nkeys <- apply(nei, 1, key)
    qidx  <- unname(idx_map[nkeys])   # 遷移先インデックス（長さ Sdim）

    # ΔE = E(q) - E(p)
    dE <- E[qidx] - E[p]

    # 受理確率 a
    if (kernel == "metropolis") {
      a <- pmin(1, exp(-beta * dE))
    } else {
      a <- 1 / (1 + exp(beta * dE))   # Glauber
    }

    # 提案分布は一様：1/Sdim
    pij <- (1 / Sdim) * a            # 長さ Sdim

    ##---------------------------
    ## 時系列制約: 明らかに過去へ戻る遷移を 0 にする
    ##---------------------------
    if (!is.null(time_state)) {
      tp <- time_state[p]
      tq <- time_state[qidx]

      # 両方NAでないものだけ判定し、「tp - tol より明らかに前」を禁止
      forbid <- !is.na(tp) & !is.na(tq) & (tq < tp - tol)

      pij[forbid] <- 0
    }

    # 自己遷移確率（残った pij の合計から計算）
    sum_pij <- sum(pij)
    if (sum_pij > 1) {
      # 数値誤差で 1 をわずかに超える場合の保険
      pij     <- pij / sum_pij
      sum_pij <- 1
    }

    p_self <- 1 - sum_pij
    if (p_self < 0) {
      # 数値誤差でマイナスに落ちたときの保険
      p_self <- 0
    }

    # triplet に自己遷移 + 近傍遷移を追加
    ii <- c(ii, rep.int(p, 1 + length(qidx)))
    jj <- c(jj, c(p, qidx))
    xx <- c(xx, c(p_self, pij))
  }

  Matrix::sparseMatrix(
    i     = ii,
    j     = jj,
    x     = xx,
    dims  = c(P, P),
    symmetric = FALSE
  )
}

# 汎用: 系列ごとの点数を数えて、線/点を自動で出し分け
add_line_and_point_layers <- function(p, df, x, y, group, color){
  # df: すでに major は factor
  library(dplyr)
  dfc <- df %>% group_by(.data[[group]]) %>% mutate(.n = n()) %>% ungroup()

  p +
    # 線は2点以上ある系列のみ
    geom_line(
      data = dfc %>% dplyr::filter(.n > 1),
      mapping = aes_string(x = x, y = y, group = group, color = color),
      linewidth = 0.7, alpha = 0.95
    ) +
    # 1点だけの系列は点で
    geom_point(
      data = dfc %>% dplyr::filter(.n == 1),
      mapping = aes_string(x = x, y = y, color = color),
      size = 2.2, alpha = 0.95
    )
}

# ----- C) Control + Cov: E(σ; ε) を epsilons 上で計算して折れ線 -----
# S_cont_cov: 各ベイシン（行）× スピン（列）の {−1,1} パターン行列
# h_cont_cov, J_cont_cov, g_cont_cov が {−1,1} 前提で整っていること（あなたの流儀どおり）
compute_curve <- function(S_mat, h, J, g, epsilons, major_vec, panel_title){
  nb <- nrow(S_mat)
  # E_mat: 行=length(epsilons), 列=ベイシン
  E_mat <- sapply(seq_len(nb), function(k){
    s <- as.numeric(S_mat[k, ])
    vapply(epsilons, function(eps) E_with_cov(s, h, J, g, eps), numeric(1))
  })
  # 単一ベイシン（ベクトル化される）に備えて行列化
  if (is.null(dim(E_mat))) E_mat <- matrix(E_mat, ncol = 1)

  df <- as.data.frame(E_mat)
  colnames(df) <- paste0("basin_", seq_len(ncol(df)))
  df$epsilon <- epsilons

  df_long <- tidyr::pivot_longer(df, -epsilon, names_to = "basin_id", values_to = "energy")

  # major を付与（ベイシン列に対応）— 長さチェックしてからリサイクル
  stopifnot(length(major_vec) == ncol(E_mat))
  df_long$major <- rep(factor(as.character(major_vec), levels = levels_all),
                       each = length(epsilons))
  df_long$panel <- panel_title
  df_long
}

time_points <- c(36, 48, 72, 96)

# 例：共変量強度 ε として正規化
epsilons    <- (time_points - min(time_points)) / (max(time_points) - min(time_points))

# Energy Calculation for cov model
E_with_cov <- function(s, h, J, g, eps) {
  # - sum_i h_i s_i
  term1 <- - sum(h * s)
  # - sum_{i<j} J_ij s_i s_j  （上三角だけ数える）
  U <- upper.tri(J, diag = FALSE)
  term2 <- - sum(J[U] * (tcrossprod(s)[U]))
  # - eps * sum_i g_i s_i
  term3 <- - eps * sum(g * s)
  term1 + term2 + term3
}

# Time Character to Numeric
.timevec <- list(
    "cont-36h" = 0,
    "cont-48h" = 0.2,
    "cont-72h" = 0.6,
    "cont-96h" = 1,
    "DAPT-36h" = 0,
    "DAPT-48h" = 0.2,
    "DAPT-72h" = 0.6,
    "DAPT-96h" = 1)

# Label Stratification
.labelStratify <- function(seurat.integrated, target.pattern){
    target <- grep(target.pattern, seurat.integrated@meta.data$sample)
    seurat.integrated[, target]
}

germlayer_colors <- c(brewer.pal(8, "Dark2")[c(2,3,1)], rgb(0,0,0, 0.5))
names(germlayer_colors) <- c("Endoderm", "Mesoderm", "Ectoderm", "NA")

# rm Hp-BMP2_4, Hp-Hox11_13b, Hp-FoxQ2, Hp-MSP130
markers <- c("Hp-Hnf6", "Hp-Chordin", "Hp-Hox7", "Hp-FoxJ1", "Hp-Bra", "Hp-Delta", "Hp-SoxC", "Hp-Awh", "Hp-Gcm", "Hp-Blimp1", "Hp-Endo16", "Hp-Ephrin", "Hp-IrxA", "Hp-Erg", "Hp-Ese", "Hp-FoxY", "Hp-Echn38", "Hp-Tph", "Hp-Pnlip-5", "Hp-Alx1", "Hp-Sm50")

fig2_markers <- c("Hp-SoxC", "Hp-Hnf6", "Hp-FoxJ1", "Hp-Chordin", "Hp-Ephrin", "Hp-Hox7", "Hp-Gcm", "Hp-Echn38", "Hp-FoxY", "Hp-Erg", "Hp-Ese", "Hp-Sm50", "Hp-Alx1", "Hp-Endo16", "Hp-Blimp1", "Hp-Bra", "Hp-IrxA", "Hp-Pnlip-5")

fig3_markers <- c("Hp-Gad", "Hp-Hdc-3", "Hp-Th", "Hp-Chat", "Hp-Ddc", "Hp-Tph", "Hp-Syt1-1-like", "Hp-Otp", "Hp-Ngn", "Hp-Ac/Sc", "Hp-Awh", "Hp-Delta", "Hp-SoxC", "Hp-Smad-ip", "Hp-Rx", "Hp-Hbn", "Hp-Nkx2.1", "Hp-Six3", "Hp-FoxQ2-1-like", "Hp-FoxQ2-1", "Hp-HesC", "Hp-Soxb1")

fig2_celltypes <- c("uncharacterized", "Pancreas", "Anus", "Endoderm", "Stomach_Intestine", "Intestine", "Stomach", "Skeleton", "Non_skeleton_mesoderm", "Blastocoelar_cell", "Germ_line_future", "Coelomic_pouch", "Pigment", "Aboral_ectoderm", "Oral_ectoderm", "Ciliary_band", "Neurons")

sample_colors <- c(brewer.pal(9, "Set1"), "black")

sample_names <- c('cont-24h', 'cont-36h', 'cont-48h', 'cont-72h', 'cont-96h', 'DAPT-24h', 'DAPT-36h', 'DAPT-48h', 'DAPT-72h', 'DAPT-96h')

group_names <- c('cont-24h', 'cont-36h', 'cont-48h', 'cont-72h', 'cont-96h', 'DAPT-24h', 'DAPT-36h', 'DAPT-48h', 'DAPT-72h', 'DAPT-96h')

conditions <- c("cont", "DAPT")

times <- c("24h", "36h", "48h", "72h", "96h")

seurat_clusters <- c(
    "3", "29", "12", "1", "17", "41", "33", "23", "38", "40",
    "25", "32", "36", "35", "39", "21", "30", "13", "43", "27",
    "14", "4", "42", "2", "19", "31", "9", "37", "15", "22",
    "8", "34", "11", "0", "16", "26", "28", "6", "5", "24",
    "20", "18", "7", "10")

blacklist_gene <- "Sp-Hrh2_3"

genes_ridgeplot_hpbase <- c("Hp-Pcna", "Hp-Srrm2-like", "Hp-Tpx2L1",
    "Hp-Lbr", "Hp-Map215prh-like", "Hp-Ndc80L")

genes_ridgeplot_echinobase <- c("Sp-Pcna", "Sp-Srrm2", "Sp-Tpx2L1",
    "Sp-Map215prh", "Sp-Ndc80L")

energy_limits <- c(-6.5, 6.5)

.firstup <- function(x){
  substr(x, 1, 1) <- toupper(substr(x, 1, 1))
  x
}

.tableToGMT <- function(gotable){
    # GeneSet's List
    gsc <- lapply(unique(gotable$TERM), function(g){
        target <- which(gotable$TERM == g)
        geneIds <- gotable$GENENAME[target]
        GeneSet(setName=g, geneIds)
    })
    GeneSetCollection(gsc)
}

.stratifySeurat <- function(seurat.integrated, group_names){
    seurat.each1 <- seurat.integrated[,
        which(seurat.integrated@meta.data$sample == group_names[1])]
    seurat.each2 <- seurat.integrated[,
        which(seurat.integrated@meta.data$sample == group_names[2])]
    seurat.each3 <- seurat.integrated[,
        which(seurat.integrated@meta.data$sample == group_names[3])]
    seurat.each4 <- seurat.integrated[,
        which(seurat.integrated@meta.data$sample == group_names[4])]
    seurat.each5 <- seurat.integrated[,
        which(seurat.integrated@meta.data$sample == group_names[5])]
    seurat.each6 <- seurat.integrated[,
        which(seurat.integrated@meta.data$sample == group_names[6])]
    seurat.each7 <- seurat.integrated[,
        which(seurat.integrated@meta.data$sample == group_names[7])]
    seurat.each8 <- seurat.integrated[,
        which(seurat.integrated@meta.data$sample == group_names[8])]
    seurat.each9 <- seurat.integrated[,
        which(seurat.integrated@meta.data$sample == group_names[9])]
    seurat.each10 <- seurat.integrated[,
        which(seurat.integrated@meta.data$sample == group_names[10])]
    list(seurat.each1, seurat.each2, seurat.each3, seurat.each4,
        seurat.each5, seurat.each6, seurat.each7, seurat.each8, seurat.each9, seurat.each10)
}

.stratifySeurat2 <- function(seurat.integrated, group_names){
    seurat.each1 <- seurat.integrated[,
        which(seurat.integrated@meta.data$sample == group_names[1])]
    seurat.each2 <- seurat.integrated[,
        which(seurat.integrated@meta.data$sample == group_names[2])]
    seurat.each3 <- seurat.integrated[,
        which(seurat.integrated@meta.data$sample == group_names[3])]
    seurat.each4 <- seurat.integrated[,
        which(seurat.integrated@meta.data$sample == group_names[4])]
    list(seurat.each1, seurat.each2, seurat.each3, seurat.each4)
}

.stratifySeurat3 <- function(seurat.integrated, group_names){
    seurat.each1 <- seurat.integrated[,
        which(seurat.integrated@meta.data$sample == group_names[1])]
    seurat.each2 <- seurat.integrated[,
        which(seurat.integrated@meta.data$sample == group_names[2])]
    seurat.each3 <- seurat.integrated[,
        which(seurat.integrated@meta.data$sample == group_names[3])]
    seurat.each4 <- seurat.integrated[,
        which(seurat.integrated@meta.data$sample == group_names[4])]
    seurat.each5 <- seurat.integrated[,
        which(seurat.integrated@meta.data$sample == group_names[5])]
    seurat.each6 <- seurat.integrated[,
        which(seurat.integrated@meta.data$sample == group_names[6])]
    seurat.each7 <- seurat.integrated[,
        which(seurat.integrated@meta.data$sample == group_names[7])]
    seurat.each8 <- seurat.integrated[,
        which(seurat.integrated@meta.data$sample == group_names[8])]
    list(seurat.each1, seurat.each2, seurat.each3, seurat.each4,
        seurat.each5, seurat.each6, seurat.each7, seurat.each8)
}

.panelPlotMeta <- function(seuratList, group_names, feature){
    gList <- list()
    
    for(i in seq_along(group_names)){
        # UMAP座標とメタデータを取得
        umap_coords <- Embeddings(seuratList[[i]], "umap")
        meta_value <- seuratList[[i]]@meta.data[[feature]]
        
        # データフレーム作成
        plot_df <- data.frame(
            UMAP_1 = umap_coords[,1],
            UMAP_2 = umap_coords[,2],
            value = meta_value
        )
        
        # プロット作成
        gList[[i]] <- ggplot(plot_df, aes(x = UMAP_1, y = UMAP_2, color = value)) +
            geom_point(size = 2) +
            scale_color_viridis_c(name = feature) +  # 連続値用のカラースケール
            labs(title = group_names[i]) +
            xlim(c(-15, 15)) + 
            ylim(c(-15, 15)) +
            theme_minimal() +
            theme(
                legend.position = "right",
                plot.title = element_text(hjust = 0.5, size = 14)
            )
    }
    
    names(gList) <- group_names
    gList <- gList[order(names(gList))]
    
    # Patch work
    (gList[[1]] | gList[[2]] | gList[[3]] | gList[[4]]) / 
    (gList[[5]] | gList[[6]] | gList[[7]] | gList[[8]])
}

.panelPlot <- function(seuratList, group_names, features){
    # Stratify
    gList <- list()
    length(gList) <- length(group_names)
    for(i in seq_along(group_names)){
        gList[[i]] <- FeaturePlot(seuratList[[i]], features=features,
            reduction = "umap", pt.size=2, label.size=6) + labs(title=group_names[i]) + xlim(c(-15,15)) + ylim(c(-15,15))

    }
    names(gList) <- group_names
    gList <- gList[order(names(gList))]
    # Patch work
    (gList[[1]] | gList[[2]] | gList[[3]] | gList[[4]] | gList[[5]]) / (gList[[6]] | gList[[7]] | gList[[8]] | gList[[9]] | gList[[10]])
}

.panelPlot2 <- function(seuratList, group_names, features){
    # Stratify
    gList <- list()
    length(gList) <- length(group_names)
    for(i in seq_along(group_names)){
        gList[[i]] <- FeaturePlot(seuratList[[i]], features=features,
            reduction = "umap", pt.size=2, label.size=6) +
            scale_color_viridis(option = "inferno", limits = c(0,1)) +
            labs(title=group_names[i]) +
            xlim(c(-15,15)) + ylim(c(-15,15))
    }
    names(gList) <- group_names
    gList <- gList[order(names(gList))]
    # Patch work
    (gList[[1]] | gList[[2]] | gList[[3]] | gList[[4]])
}

.BarPlot <- function(seurat.integrated){
    data <- seurat.integrated@meta.data[,
        c("seurat_clusters", "sample")]
    data <- .countCategory(data)
    data <- .sortByDiff(data)
    g <- ggplot(data, aes(x=seurat_clusters, y=no_cells, fill=seurat_clusters))
    g <- g + geom_bar(stat='identity')
    g <- g + xlab('Cell type')
    g <- g + ylab('# cells')
    g <- g + facet_wrap(~sample, ncol=1)
    g <- g + theme(strip.text.x = element_text(size = 25))
    g
}

.countCategory <- function(data){
    x <- as.character(unique(data[, "seurat_clusters"]))
    y <- as.character(unique(data[, "sample"]))
    xy <- as.matrix(expand.grid(x, y))
    out <- do.call("rbind", apply(xy, 1, function(z){
        targetx <- which(data[, "seurat_clusters"] == z[1])
        targety <- which(data[, "sample"] == z[2])
        no_cells <- length(intersect(targetx, targety))
        data.frame(z[1], z[2], no_cells)
    }))
    colnames(out) <- c("seurat_clusters", "sample", "no_cells")
    out
}

.sortByDiff <- function(data){
    cont_data <- data[data[, "sample"] %in% c('cont-24h', 'cont-36h', 'cont-48h', 'cont-72h', 'cont-96h'), ]
    DAPT_data <- data[data[, "sample"] %in% c('DAPT-24h', 'DAPT-36h', 'DAPT-48h', 'DAPT-72h', 'DAPT-96h'), ]
    sum_cont_data <- .sumByCelltype(cont_data)
    sum_DAPT_data <- .sumByCelltype(DAPT_data)
    diff <- sum_DAPT_data[,2] - sum_cont_data[,2]
    sum_DAPT_data[,1][order(diff)]
    data[, "seurat_clusters"] <- factor(data[, "seurat_clusters"], levels=sum_DAPT_data[,1][order(diff)])
    data
}

.sumByCelltype <- function(data){
    target <- unique(data[,1])
    no_cells <- do.call("rbind", lapply(target, function(x){
        sum(data[which(data[,1] == x), "no_cells"])
    }))
    data.frame(seurat_clusters=target, no_cells=no_cells)
}
