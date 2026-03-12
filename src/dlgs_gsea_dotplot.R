source("src/Functions.R")

args <- commandArgs(trailingOnly = TRUE)
infile_bp <- args[1]  # dlgs_gsea_bp.RData
infile_mf <- args[2]  # dlgs_gsea_mf.RData
infile_cc <- args[3]  # dlgs_gsea_cc.RData
outfile   <- args[4]  # dot plot PNG

alpha <- 0.35
max_pathways <- 15  # max pathways per ontology

load_fgsea <- function(path, ontology_label) {
  load(path)  # res_list, gos
  fg <- res_list$fg
  if (is.null(fg) || nrow(fg) == 0) return(NULL)
  sig <- fg[!is.na(fg$padj) & fg$padj < alpha, , drop = FALSE]
  if (nrow(sig) == 0) return(NULL)
  sig <- sig[order(sig$padj), ]
  if (nrow(sig) > max_pathways) sig <- sig[1:max_pathways, ]
  sig$ontology <- ontology_label
  sig
}

df_bp <- load_fgsea(infile_bp, "BP")
df_mf <- load_fgsea(infile_mf, "MF")
df_cc <- load_fgsea(infile_cc, "CC")

df <- do.call(rbind, list(df_bp, df_mf, df_cc))

if (is.null(df) || nrow(df) == 0) {
  # No significant pathways -> blank
  png(outfile, width = 10, height = 6, units = "in", res = 300)
  par(mar = c(0, 0, 0, 0))
  plot.new()
  dev.off()
  quit(save = "no")
}

# Truncate long pathway names
df$pathway_short <- ifelse(nchar(df$pathway) > 40,
                           paste0(substr(df$pathway, 1, 37), "..."),
                           df$pathway)

# Order by NES within ontology
df$pathway_short <- factor(df$pathway_short,
                           levels = rev(df$pathway_short[order(df$ontology, df$NES)]))

p <- ggplot(df, aes(x = NES, y = pathway_short)) +
  geom_point(aes(size = size, color = -log10(padj))) +
  scale_color_viridis_c(name = expression(-log[10](padj))) +
  scale_size_continuous(name = "Gene set size", range = c(2, 8)) +
  facet_grid(ontology ~ ., scales = "free_y", space = "free_y") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  labs(x = "Normalized Enrichment Score (NES)", y = NULL) +
  theme_classic(base_size = 12) +
  theme(
    text = element_text(size = 14),
    strip.text = element_text(size = 14, face = "bold"),
    strip.background = element_rect(fill = "grey90")
  )

# Dynamic height based on number of pathways
h <- max(4, nrow(df) * 0.35 + 2)
ggsave(outfile, p, width = 10, height = h, dpi = 300)
