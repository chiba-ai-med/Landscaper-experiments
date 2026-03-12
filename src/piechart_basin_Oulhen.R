source("src/Functions.R")

# Parameter
infile1 <- commandArgs(trailingOnly=TRUE)[1]
infile2 <- commandArgs(trailingOnly=TRUE)[2]
infile3 <- commandArgs(trailingOnly=TRUE)[3]
infile4 <- commandArgs(trailingOnly=TRUE)[4]
outfile <- commandArgs(trailingOnly=TRUE)[5]
# infile1 <- 'output/hpbase/integrated/seurat_annotated_landscaper.RData'
# infile2 <- 'plot/hpbase/integrated/Landscaper/Allstates_major_group.tsv'
# infile3 <- 'plot/hpbase/integrated/Landscaper/BIN_DATA'
# infile4 <- 'plot/hpbase/integrated/Landscaper/Basin.tsv'

# Load
load(infile1)
md <- seurat.integrated@meta.data

# BIN_DATA（空白区切り）
bin_df <- read.table(infile3, header=FALSE)
bin <- as.matrix(bin_df)
storage.mode(bin) <- "integer"
rownames(bin) <- rownames(md)

allst_df <- read.table(infile2, header=FALSE, sep="", stringsAsFactors=FALSE)
K <- ncol(bin)
pat_all <- as.matrix(allst_df[, 1:K])
storage.mode(pat_all) <- "integer"
state_id <- allst_df[, K+1]

# パターン照合（キー化してmatch）
key_all <- apply(pat_all, 1, paste, collapse=",")
key_bin <- apply(bin, 1, paste, collapse=",")
idx <- match(key_bin, key_all)

# Basin state_id（複数行想定）
basin_ids <- scan(infile4, what=integer(), quiet=TRUE)

cell_state <- data.frame(
  cell = rownames(md),
  celltype = md$celltype,
  state_id = state_id[idx],
  is_basin = state_id[idx] %in% basin_ids,
  stringsAsFactors = FALSE
)

# 集計：celltype ごとの割合
agg <- aggregate(is_basin ~ celltype, data = cell_state, FUN = mean)
colnames(agg)[2] <- "frac_basin"
agg$frac_non <- 1 - agg$frac_basin

# 円グラフ用にlong形式へ
pie_df <- do.call(rbind, lapply(seq_len(nrow(agg)), function(i) {
  data.frame(
    celltype = agg$celltype[i],
    kind = c("Basin", "NonBasin"),
    value = c(agg$frac_basin[i], agg$frac_non[i])
  )
}))

# 角度計算
pie_df <- pie_df[order(pie_df$celltype, pie_df$kind), ]
pie_df$start <- ave(pie_df$value, pie_df$celltype,
                    FUN = function(x) c(0, cumsum(x)[-length(x)]) * 2*pi)
pie_df$end <- pie_df$start + pie_df$value * 2*pi

basin_cols <- c("NonBasin"="grey80", "Basin"="green3")

celltype_order <- c(
  "Aboral_ectoderm",
  "Oral_ectoderm",
  "Ciliary_band",
  "Neurons"
)

pie_df$celltype <- factor(pie_df$celltype, levels = celltype_order)

# celltypeでフィルタ（指定したもののみ）
pie_df <- pie_df[!is.na(pie_df$celltype), ]

# 描画
p <- ggplot(pie_df) +
  geom_arc_bar(
    aes(
      x0 = 0, y0 = 0,
      r0 = 0, r = 1,
      start = start, end = end,
      fill = kind
    ),
    color = "white", linewidth = 0.2
  ) +
  facet_wrap(~ celltype, nrow = 1) +
  coord_fixed() +
  scale_fill_manual(values = basin_cols) +
  theme_void() +
  theme(
    strip.text = element_text(size = 10),
    legend.position = "top"
  )

ggsave(outfile, p, width = 8, height = 3, dpi = 300)