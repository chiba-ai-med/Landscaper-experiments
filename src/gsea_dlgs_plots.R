library(fgsea)
library(ggplot2)

# Parameter
args <- commandArgs(trailingOnly=TRUE)
infile  <- args[1]  # dlgs_gsea_bp.RData (res_list が入っている)
outfile <- args[2]  # FINISH など

# Load
load(infile)  # res_list, gos

# Settings
alpha <- 0.35

outdir <- dirname(outfile)

fg    <- res_list$fg
score <- res_list$score

if (!is.null(fg) && nrow(fg) > 0 && !is.null(score) && length(score) > 0) {

  sig <- fg[!is.na(fg$padj) & fg$padj < alpha, , drop = FALSE]

  if (nrow(sig) > 0) {
    for (i in seq_len(nrow(sig))) {
      pw <- sig$pathway[i]

      # pathway が gos に無い場合はスキップ
      if (is.null(gos[[pw]])) next

      # filename（安全な文字列にする）
      pw_safe <- gsub("[^A-Za-z0-9_.-]+", "_", pw)
      filename <- file.path(outdir, sprintf("dlgs__%s.png", pw_safe))

      png(file = filename, width = 800, height = 400, res = 150)

      p <- plotEnrichment(
        pathway = gos[[pw]],
        stats   = score
      )

      # 線を太く（fgseaの層は1層目が曲線）
      p$layers[[1]]$aes_params$size  <- 2.0
      p$layers[[1]]$geom_params$size <- 2.0

      # 体裁
      p <- p +
        labs(x = "Ranked genes (DLG score)", y = "Enrichment score") +
        theme_classic() +
        theme(
          text = element_text(size = 25),
          axis.title = element_text(size = 24),
          axis.text  = element_text(size = 20),
          plot.title = element_text(size = 22)
        )

      print(p)
      dev.off()
    }
  }
}

# Save (FINISH)
file.create(outfile)
