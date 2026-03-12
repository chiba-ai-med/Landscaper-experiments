source("src/Functions.R")

args <- commandArgs(trailingOnly = TRUE)
oulhen_dlgs_file <- args[1]  # Oulhen dlgs.tsv
hpbase_dlgs_file <- args[2]  # HPBase dlgs.tsv
annot_file       <- args[3]  # HpulGenome_v1_annot.xlsx
out_tsv          <- args[4]  # overlap table
out_scatter      <- args[5]  # scatter plot PNG

# --- Load DLGs ---
oulhen <- read.table(oulhen_dlgs_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
hpbase <- read.table(hpbase_dlgs_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# --- Load annotation and build mapping ---
annot <- read_excel(annot_file, sheet = 1)
annot <- as.data.frame(annot)

# Strategy 1: NR_genename -> LOC ID extraction
loc_ids <- regmatches(annot$NR_genename, regexpr("LOC\\d+", annot$NR_genename))
has_loc <- grepl("LOC\\d+", annot$NR_genename)

# LOC -> Hp name mapping
loc_to_hp <- setNames(annot$HPU_gene_name[has_loc], loc_ids)

# Strategy 2: SPU_gene_name (strip "Sp-" prefix) -> Hp name mapping
spu_valid <- annot$SPU_gene_name != "none" & !is.na(annot$SPU_gene_name)
spu_short <- sub("^Sp-", "", annot$SPU_gene_name[spu_valid])
spu_to_hp <- setNames(annot$HPU_gene_name[spu_valid], spu_short)

# --- Map Oulhen genes to HPBase gene names ---
oulhen$hp_gene <- NA_character_

# LOC match
loc_match <- match(oulhen$gene, names(loc_to_hp))
oulhen$hp_gene[!is.na(loc_match)] <- loc_to_hp[loc_match[!is.na(loc_match)]]

# SPU_gene_name match (for non-LOC genes not yet mapped)
unmapped <- is.na(oulhen$hp_gene)
spu_match <- match(oulhen$gene[unmapped], names(spu_to_hp))
oulhen$hp_gene[unmapped][!is.na(spu_match)] <- spu_to_hp[spu_match[!is.na(spu_match)]]

# Case-insensitive SPU fallback
still_unmapped <- is.na(oulhen$hp_gene)
spu_to_hp_lower <- setNames(annot$HPU_gene_name[spu_valid], tolower(spu_short))
spu_match2 <- match(tolower(oulhen$gene[still_unmapped]), names(spu_to_hp_lower))
oulhen$hp_gene[still_unmapped][!is.na(spu_match2)] <- spu_to_hp_lower[spu_match2[!is.na(spu_match2)]]

# --- Find overlap ---
hpbase_scores <- setNames(hpbase$DLG_score, hpbase$gene)
oulhen$hpbase_score <- hpbase_scores[oulhen$hp_gene]

overlap <- oulhen[!is.na(oulhen$hpbase_score), ]
overlap <- overlap[, c("gene", "hp_gene", "DLG_score", "hpbase_score")]
colnames(overlap) <- c("oulhen_gene", "hpbase_gene", "oulhen_DLG_score", "hpbase_DLG_score")
overlap <- overlap[order(-overlap$hpbase_DLG_score), ]

# --- Save overlap table ---
write.table(overlap, file = out_tsv, sep = "\t", row.names = FALSE, quote = FALSE)

# --- Scatter plot ---
cor_val <- cor(overlap$oulhen_DLG_score, overlap$hpbase_DLG_score, method = "spearman")
n_overlap <- nrow(overlap)

# Top/bottom genes to label
n_label <- 10
overlap_sorted <- overlap[order(-abs(overlap$hpbase_DLG_score) - abs(overlap$oulhen_DLG_score)), ]
top_genes <- head(overlap_sorted, n_label)

# Concordant extremes (both high positive or both high negative)
overlap$concordant <- overlap$oulhen_DLG_score * overlap$hpbase_DLG_score > 0
overlap$label <- ""
overlap$label[rownames(overlap) %in% rownames(top_genes)] <- overlap$hpbase_gene[rownames(overlap) %in% rownames(top_genes)]

p <- ggplot(overlap, aes(x = hpbase_DLG_score, y = oulhen_DLG_score)) +
  geom_point(aes(color = concordant), size = 1.2, alpha = 0.6) +
  scale_color_manual(values = c("TRUE" = "#D73027", "FALSE" = "#4575B4"),
                     labels = c("TRUE" = "Concordant", "FALSE" = "Discordant"),
                     name = "Direction") +
  geom_hline(yintercept = 0, linewidth = 0.3, linetype = "dashed") +
  geom_vline(xintercept = 0, linewidth = 0.3, linetype = "dashed") +
  geom_abline(slope = 1, intercept = 0, linewidth = 0.3, linetype = "dotted", color = "grey50") +
  ggrepel::geom_text_repel(
    aes(label = label), size = 3, max.overlaps = 20,
    segment.size = 0.2, min.segment.length = 0
  ) +
  annotate("text", x = Inf, y = -Inf,
           label = sprintf("rho = %.3f\nn = %d", cor_val, n_overlap),
           hjust = 1.1, vjust = -0.5, size = 5) +
  labs(x = "HPBase DLG score", y = "Oulhen DLG score") +
  theme_classic(base_size = 14) +
  theme(
    axis.line = element_line(linewidth = 0.8),
    axis.ticks = element_line(linewidth = 0.8)
  ) +
  coord_equal()

ggsave(out_scatter, p, width = 8, height = 8, dpi = 300)
