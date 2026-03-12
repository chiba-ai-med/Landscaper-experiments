source("src/Functions.R")

# Parameter
pc <- commandArgs(trailingOnly=TRUE)[1]
outfile <- commandArgs(trailingOnly=TRUE)[2]

######## Control ########
# Loading
basin_cont <- unlist(read.table(paste0('plot/cont/', pc, '/Landscaper/Basin.tsv')))
major_group_cont <- read.table(paste0('plot/cont/', pc, '/Landscaper/major_group.tsv'))
major_group_cont <- major_group_cont[basin_cont, ncol(major_group_cont)]
E_cont <- read.table(paste0('plot/cont/', pc, '/Landscaper/E.tsv'))
E_cont <- E_cont[basin_cont, 1]

# Create dataframe for cont (複数のtime_points値でコピー)
# time_pointsはFunctions.Rで定義済み: c(24, 36, 48, 72, 96)
df_cont <- data.frame()
for (i in seq_along(basin_cont)) {
  for (tp in time_points) {
    df_cont <- rbind(df_cont, data.frame(
      basin_id = basin_cont[i],  # Basin.tsvから読み込んだ実際の値
      time_points = tp,
      energy = E_cont[i],  # エネルギーは定数
      major = major_group_cont[i]
    ))
  }
}

######## DAPT ########
# Loading
basin_DAPT <- unlist(read.table(paste0('plot/DAPT/', pc, '/Landscaper/Basin.tsv')))
major_group_DAPT <- read.table(paste0('plot/DAPT/', pc, '/Landscaper/major_group.tsv'))
major_group_DAPT <- major_group_DAPT[basin_DAPT, ncol(major_group_DAPT)]
E_DAPT <- read.table(paste0('plot/DAPT/', pc, '/Landscaper/E.tsv'))
E_DAPT <- E_DAPT[basin_DAPT, 1]

# Create dataframe for DAPT (複数のtime_points値でコピー)
df_dapt <- data.frame()
for (i in seq_along(basin_DAPT)) {
  for (tp in time_points) {
    df_dapt <- rbind(df_dapt, data.frame(
      basin_id = basin_DAPT[i],  # Basin.tsvから読み込んだ実際の値
      time_points = tp,
      energy = E_DAPT[i],  # エネルギーは定数
      major = major_group_DAPT[i]
    ))
  }
}

######## Control + Covariate ########
h_cont_cov <- read.table(paste0('plot/cont_cov/', pc, '/Landscaper/h.tsv'))
J_cont_cov <- read.table(paste0('plot/cont_cov/', pc, '/Landscaper/J.tsv'))
g_cont_cov <- read.table(paste0('plot/cont_cov/', pc, '/Landscaper/g.txt'))
basin_cont_cov <- unlist(read.table(paste0('plot/cont_cov/', pc, '/Landscaper/Basin.tsv')))
major_group_cont_cov <- read.table(paste0('plot/cont_cov/', pc, '/Landscaper/Allstates_major_group.tsv'))
end_cont_cov <- -ncol(major_group_cont_cov)
start_cont_cov <- end_cont_cov + 1
basin_patterns_cont_cov <- major_group_cont_cov[basin_cont_cov, start_cont_cov:end_cont_cov]
major_group_cont_cov <- major_group_cont_cov[basin_cont_cov, ncol(major_group_cont_cov)]

# Preprocessing
h_cont_cov <- as.numeric(h_cont_cov[[1]])
J_cont_cov <- as.matrix(J_cont_cov)
mode(J_cont_cov) <- "numeric"
g_cont_cov <- as.numeric(g_cont_cov[[1]])
S_cont_cov <- as.matrix(basin_patterns_cont_cov)
mode(S_cont_cov) <- "numeric"

# Create dataframe for Control + Covariate
# time_pointsを正規化してepsilonとして使用 (24-96を0-1に変換)
df_cont_cov <- data.frame()
for (i in seq_along(basin_cont_cov)) {
  for (tp in time_points) {
    # time_pointsを0-1に正規化
    epsilon <- (tp - min(time_points)) / (max(time_points) - min(time_points))
    energy <- -0.5 * t(S_cont_cov[i,]) %*% J_cont_cov %*% S_cont_cov[i,] -
              sum(h_cont_cov * S_cont_cov[i,]) -
              epsilon * sum(g_cont_cov * S_cont_cov[i,])
    df_cont_cov <- rbind(df_cont_cov, data.frame(
      basin_id = basin_cont_cov[i],  # Basin.tsvから読み込んだ実際の値
      time_points = tp,
      energy = as.numeric(energy),
      major = major_group_cont_cov[i]
    ))
  }
}

######## DAPT + Covariate ########
h_DAPT_cov <- read.table(paste0('plot/DAPT_cov/', pc, '/Landscaper/h.tsv'))
J_DAPT_cov <- read.table(paste0('plot/DAPT_cov/', pc, '/Landscaper/J.tsv'))
g_DAPT_cov <- read.table(paste0('plot/DAPT_cov/', pc, '/Landscaper/g.txt'))
basin_DAPT_cov <- unlist(read.table(paste0('plot/DAPT_cov/', pc, '/Landscaper/Basin.tsv')))
major_group_DAPT_cov <- read.table(paste0('plot/DAPT_cov/', pc, '/Landscaper/Allstates_major_group.tsv'))
end_DAPT_cov <- -ncol(major_group_DAPT_cov)
start_DAPT_cov <- end_DAPT_cov + 1
basin_patterns_DAPT_cov <- major_group_DAPT_cov[basin_DAPT_cov, start_DAPT_cov:end_DAPT_cov]
major_group_DAPT_cov <- major_group_DAPT_cov[basin_DAPT_cov, ncol(major_group_DAPT_cov)]

# Preprocessing
h_DAPT_cov <- as.numeric(h_DAPT_cov[[1]])
J_DAPT_cov <- as.matrix(J_DAPT_cov)
mode(J_DAPT_cov) <- "numeric"
g_DAPT_cov <- as.numeric(g_DAPT_cov[[1]])
S_DAPT_cov <- as.matrix(basin_patterns_DAPT_cov)
mode(S_DAPT_cov) <- "numeric"

# Create dataframe for DAPT + Covariate
# time_pointsを正規化してepsilonとして使用 (24-96を0-1に変換)
df_dapt_cov <- data.frame()
for (i in seq_along(basin_DAPT_cov)) {
  for (tp in time_points) {
    # time_pointsを0-1に正規化
    epsilon <- (tp - min(time_points)) / (max(time_points) - min(time_points))
    energy <- -0.5 * t(S_DAPT_cov[i,]) %*% J_DAPT_cov %*% S_DAPT_cov[i,] -
              sum(h_DAPT_cov * S_DAPT_cov[i,]) -
              epsilon * sum(g_DAPT_cov * S_DAPT_cov[i,])
    df_dapt_cov <- rbind(df_dapt_cov, data.frame(
      basin_id = basin_DAPT_cov[i],  # Basin.tsvから読み込んだ実際の値
      time_points = tp,
      energy = as.numeric(energy),
      major = major_group_DAPT_cov[i]
    ))
  }
}

# データフレームに条件を追加して統合
df_cont$condition <- "Control"
df_dapt$condition <- "DAPT"
df_cont_cov$condition <- "Control + Covariate"
df_dapt_cov$condition <- "DAPT + Covariate"

# 全データフレームを結合
df_all <- rbind(df_cont, df_dapt, df_cont_cov, df_dapt_cov)

# conditionをfactorにして順序を設定
df_all$condition <- factor(df_all$condition,
                           levels = c("Control", "DAPT", "Control + Covariate", "DAPT + Covariate"))

# basin_idとconditionを組み合わせてユニークなグループIDを作成
df_all$group_id <- paste(df_all$condition, df_all$basin_id, sep = "_")

# legend用の表示名を作成: "Pattern X" または "Pattern X (major_group)"
# conditionとbasin_idの組み合わせごとにユニークなIDを作成（グループ分けに使用）
df_all$unique_pattern_id <- paste(df_all$condition, df_all$basin_id, sep = "_")

# 細胞型名から".数字"を除去
df_all$major_clean <- gsub("\\.\\d+$", "", df_all$major)

# 細胞型名の最初の文字を大文字に変換
df_all$major_clean <- gsub("^([a-z])", "\\U\\1", df_all$major_clean, perl = TRUE)

# 表示用のラベルはbasin_idのみを使用
df_all$legend_label <- ifelse(
  grepl("^Unknown", df_all$major_clean),  # majorが"Unknown"で始まる場合
  paste0("Pattern ", df_all$basin_id),
  paste0("Pattern ", df_all$basin_id, " (", df_all$major_clean, ")")
)

# 共有: カラーパレット（全候補を事前に割当）
# basin_idベースでユニークな組み合わせを取得（重複を除去）
legend_mapping <- df_all %>%
  select(basin_id, major_clean, legend_label) %>%
  distinct() %>%
  group_by(basin_id) %>%
  slice(1) %>%  # 各basin_idで最初の1つだけ保持
  ungroup() %>%
  arrange(basin_id)

# basin_idベースでカラーパレットを作成（同じbasin_idは同じ色）
unique_basins <- legend_mapping$basin_id
pal_basin <- setNames(scales::hue_pal()(length(unique_basins)), unique_basins)

# unique_pattern_id用のカラーパレットとラベルマップを作成
# 各unique_pattern_idに対応するbasin_idの色を割り当て
pattern_mapping <- df_all %>%
  select(unique_pattern_id, basin_id, legend_label) %>%
  distinct()

pal <- setNames(pal_basin[as.character(pattern_mapping$basin_id)], pattern_mapping$unique_pattern_id)
label_map <- setNames(pattern_mapping$legend_label, pattern_mapping$unique_pattern_id)

# legend用に重複を除去したマッピングを作成
legend_only_mapping <- df_all %>%
  select(basin_id, legend_label) %>%
  distinct() %>%
  group_by(basin_id) %>%
  slice(1) %>%
  ungroup() %>%
  arrange(basin_id)

# 4パネル共通のEnergy範囲を計算
energy_range <- range(df_all$energy, na.rm = TRUE)
energy_margin <- diff(energy_range) * 0.05
y_limits <- c(energy_range[1] - energy_margin, energy_range[2] + energy_margin)

# 共通テーマ
common_theme <- theme_classic(base_size = 20) +
  theme(legend.position = "none")

# 統合データフレームから各条件のプロットを作成
p_cont <- ggplot(df_all[df_all$condition == "Control",],
                 aes(x = time_points, y = energy, color = unique_pattern_id, group = group_id)) +
  geom_line(linewidth = 2.8, alpha = 0.95) +
  geom_point(size = 5, alpha = 0.95) +
  scale_color_manual(values = pal, labels = label_map, drop = FALSE) +
  scale_y_continuous(limits = y_limits) +
  labs(title = "Control", x = "Time (Hour)", y = "Energy") +
  common_theme

p_dapt <- ggplot(df_all[df_all$condition == "DAPT",],
                 aes(x = time_points, y = energy, color = unique_pattern_id, group = group_id)) +
  geom_line(linewidth = 2.8, alpha = 0.95) +
  geom_point(size = 5, alpha = 0.95) +
  scale_color_manual(values = pal, labels = label_map, drop = FALSE) +
  scale_y_continuous(limits = y_limits) +
  labs(title = "DAPT", x = "Time (Hour)", y = "Energy") +
  common_theme

p_cont_cov <- ggplot(df_all[df_all$condition == "Control + Covariate",],
                     aes(x = time_points, y = energy, color = unique_pattern_id, group = group_id)) +
  geom_line(linewidth = 2.8, alpha = 0.95) +
  geom_point(size = 5, alpha = 0.95) +
  scale_color_manual(values = pal, labels = label_map, drop = FALSE) +
  scale_y_continuous(limits = y_limits) +
  labs(title = "Control + Covariate", x = "Time (Hour)", y = "Energy") +
  common_theme

p_dapt_cov <- ggplot(df_all[df_all$condition == "DAPT + Covariate",],
                     aes(x = time_points, y = energy, color = unique_pattern_id, group = group_id)) +
  geom_line(linewidth = 2.8, alpha = 0.95) +
  geom_point(size = 5, alpha = 0.95) +
  scale_color_manual(values = pal, labels = label_map, drop = FALSE) +
  scale_y_continuous(limits = y_limits) +
  labs(title = "DAPT + Covariate", x = "Time (Hour)", y = "Energy") +
  common_theme

# ---- 凡例用のダミープロットを作成 ----
# 重複を除去した凡例を作成（basin_idごとに1つ）
# basin_idベースでダミーデータを作成
dummy_data <- data.frame(
  x = 0.5,
  y = 0,
  basin_id = legend_only_mapping$basin_id,
  label = legend_only_mapping$legend_label
)

p_legend <- ggplot(dummy_data, aes(x = x, y = y, color = as.factor(basin_id), group = basin_id)) +
  geom_line(linewidth = 2.8, alpha = 0.95) +
  scale_color_manual(values = pal_basin, labels = setNames(dummy_data$label, dummy_data$basin_id), drop = FALSE) +
  theme_classic(base_size = 20) +
  theme(legend.position = "right",
        legend.title = element_blank())

# 凡例だけを抽出
library(cowplot)
legend <- get_legend(p_legend)

# ---- 2×2 に結合（凡例を追加） ----
combined_plots <- (p_cont + p_dapt) / (p_cont_cov + p_dapt_cov)
combined <- plot_grid(combined_plots, legend, rel_widths = c(1, 0.3))

ggsave(outfile, combined, width = 15, height = 9, dpi = 300)

