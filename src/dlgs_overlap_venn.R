source("src/Functions.R")

args <- commandArgs(trailingOnly = TRUE)
overlap_tsv <- args[1]  # dlgs_overlap_{comparison}.tsv
oulhen_dlgs <- args[2]  # Oulhen dlgs.tsv
hpbase_dlgs <- args[3]  # HPBase dlgs.tsv
out_top     <- args[4]  # Venn diagram PNG (positive DLGs)
out_bottom  <- args[5]  # Venn diagram PNG (negative DLGs)

pct <- 0.10  # top/bottom 10%

# --- Load data ---
overlap <- read.table(overlap_tsv, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
oulhen  <- read.table(oulhen_dlgs, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
hpbase  <- read.table(hpbase_dlgs, header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# Use only mappable genes for fair comparison
mappable_hp <- unique(overlap$hpbase_gene)

# Oulhen: map to hpbase gene names via overlap table
oulhen_map <- setNames(overlap$hpbase_gene, overlap$oulhen_gene)
oulhen$hp_gene <- oulhen_map[oulhen$gene]
oulhen_mapped <- oulhen[!is.na(oulhen$hp_gene), ]
oulhen_mapped <- oulhen_mapped[order(-oulhen_mapped$DLG_score), ]

hpbase_mapped <- hpbase[hpbase$gene %in% mappable_hp, ]
hpbase_mapped <- hpbase_mapped[order(-hpbase_mapped$DLG_score), ]

# --- Percentile-based selection (top/bottom 10% of each dataset) ---
n_oulhen_top <- ceiling(nrow(oulhen_mapped) * pct)
n_hpbase_top <- ceiling(nrow(hpbase_mapped) * pct)

oulhen_top <- head(oulhen_mapped$hp_gene, n_oulhen_top)
hpbase_top <- head(hpbase_mapped$gene, n_hpbase_top)

oulhen_bot <- tail(oulhen_mapped$hp_gene, n_oulhen_top)
hpbase_bot <- tail(hpbase_mapped$gene, n_hpbase_top)

# --- Venn diagram drawing function (grid) ---
draw_venn <- function(set_a, set_b, label_a, label_b, title, outfile) {
  only_a <- length(setdiff(set_a, set_b))
  only_b <- length(setdiff(set_b, set_a))
  both   <- length(intersect(set_a, set_b))
  shared_genes <- sort(intersect(set_a, set_b))

  png(outfile, width = 8, height = 8, units = "in", res = 300)

  grid.newpage()

  # Title
  grid.text(title, x = 0.5, y = 0.95, gp = gpar(fontsize = 18, fontface = "bold"))

  # Subtitle with set sizes
  grid.text(sprintf("Oulhen: %d genes, HPBase: %d genes", length(set_a), length(set_b)),
            x = 0.5, y = 0.90, gp = gpar(fontsize = 12, col = "grey40"))

  # Circles
  grid.circle(x = 0.35, y = 0.52, r = 0.22,
              gp = gpar(fill = alpha("#D73027", 0.3), col = "#D73027", lwd = 2))
  grid.circle(x = 0.65, y = 0.52, r = 0.22,
              gp = gpar(fill = alpha("#4575B4", 0.3), col = "#4575B4", lwd = 2))

  # Labels
  grid.text(label_a, x = 0.18, y = 0.80, gp = gpar(fontsize = 14, fontface = "bold", col = "#D73027"))
  grid.text(label_b, x = 0.82, y = 0.80, gp = gpar(fontsize = 14, fontface = "bold", col = "#4575B4"))

  # Counts
  grid.text(only_a, x = 0.25, y = 0.52, gp = gpar(fontsize = 20, fontface = "bold"))
  grid.text(both,   x = 0.50, y = 0.52, gp = gpar(fontsize = 20, fontface = "bold"))
  grid.text(only_b, x = 0.75, y = 0.52, gp = gpar(fontsize = 20, fontface = "bold"))

  # Shared gene names (bottom)
  if (length(shared_genes) > 0) {
    gene_text <- paste(shared_genes, collapse = ", ")
    wrapped <- strwrap(gene_text, width = 80)
    max_lines <- min(length(wrapped), 8)
    wrapped <- wrapped[seq_len(max_lines)]
    if (length(strwrap(gene_text, width = 80)) > max_lines) {
      wrapped[max_lines] <- paste0(wrapped[max_lines], " ...")
    }
    for (i in seq_along(wrapped)) {
      grid.text(wrapped[i], x = 0.5, y = 0.22 - (i - 1) * 0.03,
                gp = gpar(fontsize = 8))
    }
    grid.text("Shared genes:", x = 0.5, y = 0.26, gp = gpar(fontsize = 10, fontface = "bold"))
  }

  dev.off()
}

# --- Draw ---
draw_venn(oulhen_top, hpbase_top,
          sprintf("Oulhen (top %d%%)", pct * 100),
          sprintf("HPBase (top %d%%)", pct * 100),
          sprintf("Top %d%% DLGs overlap", pct * 100),
          out_top)

draw_venn(oulhen_bot, hpbase_bot,
          sprintf("Oulhen (bottom %d%%)", pct * 100),
          sprintf("HPBase (bottom %d%%)", pct * 100),
          sprintf("Bottom %d%% DLGs overlap", pct * 100),
          out_bottom)
