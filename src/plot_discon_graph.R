# plot_discon_graph.R
suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
})

source("src/Functions.R")

args <- commandArgs(trailingOnly = TRUE)
infile   <- args[1]
outfile1 <- args[2]  # 直線スタイル
outfile2 <- args[3]  # 直角スタイル

# Load objects: hc, dend, node_order, df_xy, subtrees, leaves, paths, dg_skeleton
load(infile)

if (is.null(hc) || is.null(dg_skeleton) || is.null(df_xy)) {
  # Empty (単一ベイシンなど)
  file.create(outfile1)
  file.create(outfile2)
  quit(save="no")
}

# ------------------------------------------------------------
# ラベル用データ
# ・葉（basin）のラベル：dg_skeleton$state（既存）
# ・分岐点（transition）のラベル：df_xy$tipping_label（あれば）
# ------------------------------------------------------------

# Basinラベル（xend,yend は葉の座標）
lab_basin <- dg_skeleton %>%
  filter(!is.na(state)) %>%
  transmute(x = xend, y = yend - 0.05, label = as.character(state))

# Transitionラベル（tipping_label があれば重ねる）
lab_tip <- df_xy %>%
  filter(state == "transition", !is.na(tipping_label), tipping_label != "") %>%
  transmute(x = x, y = y + 0.05, label = tipping_label)

# ------------------------------------------------------------
# 1) 直線スタイル
# ------------------------------------------------------------
g1 <- ggplot(dg_skeleton) +
  geom_segment(aes(x = x, y = y, xend = xend, yend = yend)) +
  # ラベル：葉
  geom_text(data = lab_basin, aes(x = x, y = y, label = label), size = 3) +
  # ラベル：分岐点（tipping）
  geom_text(data = lab_tip, aes(x = x, y = y, label = label), size = 3) +
  theme_bw() +
  theme(
    axis.text.x  = element_blank(),
    axis.ticks.x = element_blank(),
    panel.grid   = element_blank()
  ) +
  labs(x = "", y = "Energy")

# ------------------------------------------------------------
# 2) 直角スタイル
# ------------------------------------------------------------
g2 <- ggplot(dg_skeleton) +
  geom_segment(aes(x = x,    y = y, xend = xend, yend = y)) +
  geom_segment(aes(x = xend, y = y, xend = xend, yend = yend)) +
  # ラベル：葉
  geom_text(data = lab_basin, aes(x = x, y = y, label = label), size = 3) +
  # ラベル：分岐点（tipping）
  geom_text(data = lab_tip, aes(x = x, y = y, label = label), size = 3) +
  theme_bw() +
  theme(
    axis.text.x  = element_blank(),
    axis.ticks.x = element_blank(),
    panel.grid   = element_blank()
  ) +
  labs(x = "", y = "Energy")

# ------------------------------------------------------------
# Save
# ------------------------------------------------------------
ggsave(file = outfile1, plot = g1, height = 6, width = 8)
ggsave(file = outfile2, plot = g2, height = 6, width = 8)
