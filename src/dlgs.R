source("src/Functions.R")

args <- commandArgs(trailingOnly = TRUE)

# Input
cont_seurat_file <- args[1]   # output/cont/{pc}/seurat_landscaper.RData
dapt_seurat_file <- args[2]   # output/DAPT/{pc}/seurat_landscaper.RData

# Output
output_file <- args[3]        # output/cont_DAPT/dlgs/{pc}/dlgs.tsv

# Dimensionality (PCs)
K <- as.integer(args[4])

##----------------------------------------------------------------------
## Load
##----------------------------------------------------------------------
load(cont_seurat_file)   # loads seurat.integrated
seurat.cont <- seurat.integrated

load(dapt_seurat_file)   # loads seurat.integrated
seurat.dapt <- seurat.integrated

##----------------------------------------------------------------------
## BIN_DATA (cell x PC, {-1,+1})
##----------------------------------------------------------------------
BIN_cont <- seurat.cont@misc$landscaper$BIN_DATA
BIN_dapt <- seurat.dapt@misc$landscaper$BIN_DATA

##----------------------------------------------------------------------
## Delta_m (PC-wise mean shift in binary states)
##----------------------------------------------------------------------
m_cont <- colMeans(BIN_cont)
m_dapt <- colMeans(BIN_dapt)
delta_m <- as.numeric(m_dapt - m_cont)

## Normalize (direction only)
delta_m_unit <- delta_m / sqrt(sum(delta_m^2))

##----------------------------------------------------------------------
## Gene loadings (gene x PC) from Seurat PCA
##----------------------------------------------------------------------
W <- Seurat::Loadings(seurat.cont, reduction = "pca")
W <- W[, 1:K, drop = FALSE]

##----------------------------------------------------------------------
## DLG score: cosine( w_g , delta_m )
##----------------------------------------------------------------------
wg_norm <- sqrt(rowSums(W^2))
score <- as.numeric((W %*% delta_m_unit)) / wg_norm

dlgs_df <- data.frame(
    gene = rownames(W),
    DLG_score = score,
    stringsAsFactors = FALSE
)

dlgs_df <- dlgs_df[order(dlgs_df$DLG_score, decreasing = TRUE), ]

##----------------------------------------------------------------------
## Save
##----------------------------------------------------------------------
write.table(
    dlgs_df,
    file = output_file,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    col.names = TRUE
)
