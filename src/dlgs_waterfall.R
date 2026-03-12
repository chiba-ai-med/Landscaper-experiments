source("src/Functions.R")

args <- commandArgs(trailingOnly = TRUE)
infile  <- args[1]  # dlgs.tsv
outfile <- args[2]  # waterfall plot PNG

n_label <- 15  # number of top/bottom genes to label

dlgs <- read.table(infile, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
dlgs <- dlgs[order(dlgs$DLG_score, decreasing = TRUE), ]
dlgs$rank <- seq_len(nrow(dlgs))

# Separate top and bottom gene labels
n_top <- min(n_label, nrow(dlgs))
n_genes <- nrow(dlgs)

dlgs$label_top <- ""
dlgs$label_top[seq_len(n_top)] <- dlgs$gene[seq_len(n_top)]

dlgs$label_bot <- ""
dlgs$label_bot[seq(n_genes - n_top + 1, n_genes)] <- dlgs$gene[seq(n_genes - n_top + 1, n_genes)]

# Color by sign
dlgs$direction <- ifelse(dlgs$DLG_score > 0, "positive", "negative")

y_range <- range(dlgs$DLG_score)
y_expand <- diff(y_range) * 0.4

p <- ggplot(dlgs, aes(x = rank, y = DLG_score, fill = direction)) +
  geom_bar(stat = "identity", width = 1) +
  scale_fill_manual(values = c("positive" = "#D73027", "negative" = "#4575B4"),
                    guide = "none") +
  geom_hline(yintercept = 0, linewidth = 0.3) +
  ggrepel::geom_text_repel(
    aes(label = label_top),
    size = 3.6, max.overlaps = 30, segment.size = 0.2,
    min.segment.length = 0,
    nudge_x = n_genes * 0.25,
    direction = "y",
    xlim = c(1, n_genes),
    ylim = c(y_range[1] - y_expand, y_range[2] + y_expand),
    force = 2, force_pull = 0.3
  ) +
  ggrepel::geom_text_repel(
    aes(label = label_bot),
    size = 3.6, max.overlaps = 30, segment.size = 0.2,
    min.segment.length = 0,
    nudge_x = -n_genes * 0.25,
    direction = "y",
    ylim = c(y_range[1] - y_expand, y_range[2] + y_expand),
    xlim = c(1, n_genes),
    force = 2, force_pull = 0.3
  ) +
  scale_x_continuous(limits = c(1, n_genes), expand = c(0, 0)) +
  scale_y_continuous(limits = c(y_range[1] - y_expand, y_range[2] + y_expand), expand = c(0, 0)) +
  labs(x = "Gene rank", y = "DLG score") +
  theme_classic(base_size = 14) +
  theme(
    text = element_text(size = 22),
    axis.line = element_line(linewidth = 1.0),
    axis.ticks = element_line(linewidth = 1.0),
    plot.margin = margin(10, 10, 10, 10, unit = "mm")
  )

ggsave(outfile, p, width = 8, height = 8, dpi = 300)
