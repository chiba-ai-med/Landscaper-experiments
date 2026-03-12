library(Seurat)
library(SingleCellExperiment)
library(schex)
library(ggplot2)

# Parameter
args <- commandArgs(trailingOnly=TRUE)
infile1 <- args[1]  # basin_gsea_{go}.RData
infile2 <- args[2]  # seurat_landscaper.RData
outfile <- args[3]  # FINISH

# Settings
alpha <- 0.35
nbins <- 15
max_genes_per_pathway <- 30

# Load GSEA results
load(infile1)  # res_list, gos

# Load Seurat object
load(infile2)  # seurat.integrated

# Seurat v5 compatibility
if (inherits(seurat.integrated[["RNA"]], "Assay5")) {
  seurat.integrated[["RNA"]] <- JoinLayers(seurat.integrated[["RNA"]])
}

# Convert to SingleCellExperiment and create hexbin
sce <- as.SingleCellExperiment(seurat.integrated)
sce <- make_hexbin(sce, nbins = nbins, dimension_reduction = "UMAP")

# Available genes
available_genes <- rownames(sce)

outdir <- dirname(outfile)

for (b in seq_along(res_list)) {
  fg    <- res_list[[b]]$fg
  score <- res_list[[b]]$score

  if (is.null(fg) || nrow(fg) == 0) next
  if (is.null(score) || length(score) == 0) next

  sig <- fg[!is.na(fg$padj) & fg$padj < alpha, , drop = FALSE]
  if (nrow(sig) == 0) next

  for (i in seq_len(nrow(sig))) {
    pw <- sig$pathway[i]
    if (is.null(gos[[pw]])) next

    # Leading edge genes
    le_genes <- sig$leadingEdge[[i]]
    if (is.null(le_genes) || length(le_genes) == 0) next

    # Filter to genes present in expression matrix
    le_genes <- intersect(le_genes, available_genes)
    if (length(le_genes) == 0) next

    # Cap at max_genes_per_pathway
    if (length(le_genes) > max_genes_per_pathway) {
      le_genes <- le_genes[seq_len(max_genes_per_pathway)]
    }

    pw_safe <- gsub("[^A-Za-z0-9_.-]+", "_", pw)

    for (gene in le_genes) {
      gene_safe <- gsub("[^A-Za-z0-9_.-]+", "_", gene)
      filename <- file.path(
        outdir,
        sprintf("basin_%02d__%s__%s.png", b, pw_safe, gene_safe)
      )

      png(file = filename, width = 700, height = 700)
      p <- plot_hexbin_feature(
        sce,
        action  = "median",
        type    = "logcounts",
        feature = gene
      ) +
        ggtitle(sprintf("Basin %d / %s / %s", b, pw, gene)) +
        theme(plot.title = element_text(size = 10))
      print(p)
      dev.off()
    }
  }
}

# FINISH marker
file.create(outfile)
