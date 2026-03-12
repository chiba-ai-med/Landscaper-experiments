library(fgsea)
library(dplyr)
library(GSEABase)

# Parameter
args <- commandArgs(trailingOnly=TRUE)
infile1 <- args[1]  # dlgs.tsv
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

# Load DLGs
dlgs <- read.table(infile1, header=TRUE, sep="\t", stringsAsFactors=FALSE)

# Create named vector for fgsea (gene -> DLG_score)
score <- setNames(dlgs$DLG_score, dlgs$gene)
score <- score[is.finite(score)]
score <- sort(score, decreasing = TRUE)

# Load GO
objname <- load(infile2)
gmt_go <- get(objname)
rm(list = objname)
gos <- to_fgsea_pathways(gmt_go)

# Run fgsea
set.seed(1)
fg <- fgsea(
  pathways = gos,
  stats = score,
  minSize = 5,
  maxSize = 500,
  nperm = 10000
) %>% arrange(padj, desc(abs(NES)))

# Save
res_list <- list(
  score = score,
  fg = fg
)

save(res_list, gos, file = outfile)
