library(Seurat)
library(fgsea)
library(dplyr)
library(Matrix)
library(GSEABase)

# Parameter
args <- commandArgs(trailingOnly=TRUE)
infile1 <- args[1]  # seurat_landscaper.RData (for expression, state_id, Basin)
infile2 <- args[2]  # go_*.RData
outfile <- args[3]

# Gene sets -> fgsea pathways
to_fgsea_pathways <- function(gmt_go) {
  gs_list <- if (inherits(gmt_go, "GeneSetCollection")) as.list(gmt_go) else gmt_go
  pathways <- lapply(gs_list, geneIds)
  nm <- vapply(gs_list, setName, character(1))
  nm[nm == ""] <- vapply(gs_list, setIdentifier, character(1))[nm == ""]
  names(pathways) <- nm
  pathways
}

# Load Seurat
load(infile1)
cell_state_id <- seurat.integrated$state_id
Basin <- seurat.integrated@misc$landscaper$Basin
basin_ids <- matrix(Basin, ncol = 1)
storage.mode(basin_ids) <- "integer"

# Load GO
objname <- load(infile2)
gmt_go <- get(objname)
rm(list = objname)
gos <- to_fgsea_pathways(gmt_go)

# Expression matrix (join layers for Seurat v5 compatibility)
if (inherits(seurat.integrated[["RNA"]], "Assay5")) {
  seurat.integrated[["RNA"]] <- JoinLayers(seurat.integrated[["RNA"]])
}
expr <- GetAssayData(seurat.integrated, assay="RNA", layer="data")

# ---- BasinごとにGSEA ----
res_list <- vector("list", nrow(basin_ids))

for (b in seq_len(nrow(basin_ids))) {
  ids <- basin_ids[b,]
  ids <- ids[is.finite(ids)]
  ids <- ids[ids > 0]

  is_basin <- cell_state_id %in% ids

  mu_basin    <- Matrix::rowMeans(expr[, is_basin, drop=FALSE])
  mu_nonbasin <- Matrix::rowMeans(expr[, !is_basin, drop=FALSE])

  score <- mu_basin - mu_nonbasin
  score <- score[is.finite(score)]
  score <- sort(score, decreasing = TRUE)

  set.seed(1)
  fg <- fgsea(
    pathways = gos,
    stats = score,
    minSize = 5,
    maxSize = 500,
    nperm = 10000
  ) %>% arrange(padj, desc(abs(NES)))

  res_list[[b]] <- list(
    basin_ids = ids,
    n_basin = sum(is_basin),
    n_nonbasin = sum(!is_basin),
    score = score,
    fg = fg
  )
}

names(res_list) <- paste0("basin_", seq_along(res_list))

save(res_list, gos, file = outfile)
