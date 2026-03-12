source("src/Functions.R")

# Parameter
infile1 <- commandArgs(trailingOnly=TRUE)[1]
infile2 <- commandArgs(trailingOnly=TRUE)[2]
infile3 <- commandArgs(trailingOnly=TRUE)[3]
infile4 <- commandArgs(trailingOnly=TRUE)[4]
outfile <- commandArgs(trailingOnly=TRUE)[5]
# infile1 <- 'output/hpbase/DAPT_cov/seurat_annotated_landscaper.RData'
# infile2 <- 'plot/hpbase/DAPT_cov/Landscaper/Allstates_major_group.tsv'
# infile3 <- 'plot/hpbase/DAPT_cov/Landscaper/BIN_DATA'
# infile4 <- 'plot/hpbase/DAPT_cov/Landscaper/Basin.tsv'

# -----------------------------
# Load Seurat
# -----------------------------
load(infile1)
if (!exists("seurat.integrated")) stop("RData内に seurat.integrated がありません")
md <- seurat.integrated@meta.data

# time抽出: cont-36h -> 36
md$time <- as.integer(sub(".*-([0-9]+)h$", "\\1", md$sample))
if (any(is.na(md$time))) stop("time抽出に失敗: md$sample の形式を確認して下さい")

# -----------------------------
# Read BIN_DATA (space separated)
# -----------------------------
bin_df <- read.table(infile3, header = FALSE, sep = "", stringsAsFactors = FALSE)
bin <- as.matrix(bin_df)
storage.mode(bin) <- "integer"

if (nrow(bin) != nrow(md)) {
  stop(sprintf("BIN_DATA行数(%d) != Seurat meta行数(%d)。行順対応の前提が崩れています。",
               nrow(bin), nrow(md)))
}
rownames(bin) <- rownames(md)

# -----------------------------
# Read Allstates_major_group.tsv
# -----------------------------
allst_df <- read.table(infile2, header = FALSE, sep = "", stringsAsFactors = FALSE)
K <- ncol(bin)
if (ncol(allst_df) < K + 2) stop("Allstates_major_group.tsv の列数が足りません（K+2未満）")

pat_all <- as.matrix(allst_df[, 1:K])
storage.mode(pat_all) <- "integer"
state_id <- allst_df[, K + 1]
# major_group <- allst_df[, K + 2]  # ここでは使わない

# パターン照合（キー化してmatch）
key_all <- apply(pat_all, 1, paste, collapse = ",")
key_bin <- apply(bin, 1, paste, collapse = ",")
idx <- match(key_bin, key_all)
if (anyNA(idx)) {
  bad <- sum(is.na(idx))
  stop(sprintf("BIN_DATAの%d行がAllstatesにマッチしません（K不一致/値の不一致の可能性）。", bad))
}

# Basin state_id（複数行想定）
basin_ids <- scan(infile4, what = integer(), quiet = TRUE)
if (length(basin_ids) == 0) stop("Basin.tsv が空です")

cell_state <- data.frame(
  cell = rownames(md),
  time = md$time,
  celltype = md$celltype,
  state_id = state_id[idx],
  is_basin = state_id[idx] %in% basin_ids,
  stringsAsFactors = FALSE
)

# -----------------------------
# glm を RData に保存
# -----------------------------
# 1) シンプル：全体で time の効果（celltypeを共変量に入れる）
#   ※ celltype横断的に「時間でBasin割合が変わるか」を一括で見る
cell_state$time <- as.numeric(cell_state$time)
cell_state$celltype <- factor(cell_state$celltype)

fit_glm <- glm(is_basin ~ time + celltype, data = cell_state, family = binomial())
fit_glm_summary <- summary(fit_glm)

# 2) 追加で「celltypeごとに傾きが違うか」を見たければ interaction
fit_glm_int <- glm(is_basin ~ time * celltype, data = cell_state, family = binomial())
anova_time_by_celltype <- anova(fit_glm, fit_glm_int, test = "Chisq")

save(
  cell_state,
  fit_glm,
  fit_glm_summary,
  fit_glm_int,
  anova_time_by_celltype,
  file = outfile)
