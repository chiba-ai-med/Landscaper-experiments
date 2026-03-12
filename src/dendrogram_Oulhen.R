source("src/Functions.R")

args <- commandArgs(trailingOnly = TRUE)
infile1 <- args[1]
infile2 <- args[2]
infile3 <- args[3]
infile4 <- args[4]
infile5 <- args[5]
outfile <- args[6]
# infile1 <- 'plot/hpbase/integrated/Landscaper/E.tsv'
# infile2 <- 'plot/hpbase/integrated/Landscaper/Basin.tsv'
# infile3 <- 'plot/hpbase/integrated/Landscaper/EnergyBarrier.tsv'
# infile4 <- 'plot/hpbase/integrated/Landscaper/Allstates_major_group.tsv'
# infile5 <- 'plot/hpbase/integrated/Landscaper/TippingState.tsv'

# ---- Load ----
E <- unlist(read.table(infile1, header = FALSE))
Basin <- unlist(read.table(infile2, header = FALSE))
EnergyBarrier <- as.matrix(read.table(infile3, header = FALSE))
if (file.size(infile4) != 0) {
  Group <- read.table(infile4, header = FALSE, stringsAsFactors = FALSE)
} else {
  Group <- NULL
}
TippingState <- as.matrix(read.table(infile5, header = FALSE))

# ---- 単一ベイシンなら空オブジェクトで保存して終了 ----
if (length(Basin) == 1) {
  hc <- NULL; dend <- NULL; node_order <- NULL; df_xy <- NULL
  subtrees <- NULL; leaves <- NULL; paths <- NULL; dg_skeleton <- NULL
  save(hc, dend, node_order, df_xy, subtrees, leaves, paths, dg_skeleton, file = outfile)
  quit(save = "no")
}

# ---- Dendrogram 構築 ----
state_energy <- data.frame(state = seq_along(E) - 1L, E = E)

hc <- EnergyBarrier |>
  as.dist() |>
  hclust(method = "single")

dend <- as.dendrogram(hc)

# ノード順（葉の並び）
node_order <- tibble(leaf = hc$order, count = seq_along(hc$order))

# ノード座標の抽出（Functions.R の get_nodes_xy() を使用）
df_xy <- get_nodes_xy(dend) %>%
  data.frame() %>%
  as_tibble() %>%
  magrittr::set_colnames(c("x", "y")) %>%
  mutate(node  = seq_len(nrow(.)),
         state = if_else(y == 0, "basin", "transition")) %>%
  group_by(state) %>%
  mutate(count = row_number()) %>%
  left_join(node_order, by = "count") %>%
  mutate(count = if_else(state == "transition", count, leaf),
         node_name = str_c(state, count, sep = "_")) %>%
  ungroup() %>%
  # basin ノードの y に E を上書き
  left_join(
    state_energy |>
      filter(state %in% (Basin - 1L)) |>
      mutate(stateID = match(state, Basin - 1L),
             stateID = paste0("basin_", stateID)) |>
      select(-state),
    by = c("node_name" = "stateID")
  ) %>%
  mutate(y = if_else(is.na(E), y, E))

# 葉番号(=hcの順) とベイシン state の対応
leaf2state <- tibble(state = Basin) |>
  mutate(leaf_no = row_number()) |>
  mutate(leaf_no = as.character(leaf_no))
# 表示ラベルとして Group の最後の列を使う（あれば）
if (!is.null(Group)) {
  leaf2state$state <- Group[leaf2state$state, ncol(Group)]
}

# サブツリー分割 & 経路取り（Functions.R の partition_leaves(), pathRoutes() を使用）
subtrees <- partition_leaves(dend)
leaves   <- subtrees[[1]]
paths <- lapply(leaves, function(x, subtrees) pathRoutes(x, subtrees), subtrees = subtrees)

# 骨格エッジ
dg_skeleton <- tibble(
  leaf = rep(seq_along(paths), times = lengths(paths)),
  from = flatten_int(paths)
) |>
  mutate(to = lead(from, 1)) |>
  filter(to > from) |>
  select(from, to) |>
  distinct() %>%
  left_join(df_xy %>% select(x, y, node, node_name), by = c("to" = "node")) %>%
  dplyr::rename(xend = x, yend = y) %>%
  left_join(df_xy %>% select(x, y, node), by = c("from" = "node")) %>%
  mutate(leaf    = if_else(str_detect(node_name, "basin"), node_name, NA_character_),
         leaf_no = str_split_i(leaf, "_", 2)) %>%
  left_join(leaf2state, by = "leaf_no")

# ============================================================
#  分岐点（transitionノード）に tipping state 由来ラベルを付与
# ============================================================
# hc$merge から，内部ノードkの左右の葉集合を返す補助（hclust 仕様）
get_leaves_sets <- function(merge_mat, k) {
  left  <- merge_mat[k, 1]; right <- merge_mat[k, 2]
  get_leaf_set <- function(x) {
    if (x < 0) return(-x) else return(unlist(get_leaves_sets(merge_mat, x)))
  }
  list(L = get_leaf_set(left), R = get_leaf_set(right))
}

# leaf_no -> ベイシンのインデックス i（EnergyBarrier/TippingState の添字 1..n）
leaf2basin_idx <- tibble(leaf_no = as.character(seq_along(Basin)),
                         i = seq_along(Basin))

# 追加列（分岐点ラベル）
df_xy$tipping_label <- NA_character_

# 一致許容（高さ = バリアの数値誤差）
eps_h <- 1e-8
nB <- length(Basin)
internal_nodes <- seq_len(nB - 1)  # hclust の内部ノード番号
node_height <- hc$height

for (k in internal_nodes) {
  # 内部ノード k の左右葉集合（leaf index 1..nB）
  lr <- get_leaves_sets(hc$merge, k)
  L <- lr$L; R <- lr$R
  if (length(L) == 0 || length(R) == 0) next

  # それぞれをベイシン添字 i に変換
  Li <- leaf2basin_idx$i[ match(as.character(L), leaf2basin_idx$leaf_no) ]
  Ri <- leaf2basin_idx$i[ match(as.character(R), leaf2basin_idx$leaf_no) ]
  Li <- Li[!is.na(Li)]; Ri <- Ri[!is.na(Ri)]
  if (length(Li) == 0 || length(Ri) == 0) next

  # この内部ノードの高さ h
  h <- node_height[k]

  # (i,j) 候補のうち、EnergyBarrier[i,j] が h に一致（±eps）するもの
  cand <- expand.grid(i = Li, j = Ri)
  if (nrow(cand) == 0) next
  eb_vals <- mapply(function(i, j) EnergyBarrier[i, j], cand$i, cand$j)
  sel <- which(abs(eb_vals - h) <= eps_h)
  if (length(sel) == 0) sel <- which.min(abs(eb_vals - h))  # fallback: 最も近いもの

  i_sel <- cand$i[sel[1]]; j_sel <- cand$j[sel[1]]

  # tipping state（0始まり）→ +1 で行番号
  s_star0 <- TippingState[i_sel, j_sel]
  tip_lab <- NA_character_
  if (!is.na(s_star0)) {
    s_star <- s_star0 + 1L
    if (!is.null(Group)) {
      tip_lab <- as.character(Group[s_star, ncol(Group)])
    } else {
      tip_lab <- paste0("state_", s_star0)  # Groupが無ければ状態ID表示
    }
  }

  # df_xy の該当 transition ノードにラベルを入れる
  node_name_k <- paste0("transition_", k)
  df_xy$tipping_label[df_xy$node_name == node_name_k] <- tip_lab
}

# ---- Save ----
save(hc, dend, node_order, df_xy, subtrees, leaves, paths, dg_skeleton, file = outfile)
