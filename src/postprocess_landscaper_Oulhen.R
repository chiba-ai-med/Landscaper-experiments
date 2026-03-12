source("src/Functions.R")

args <- commandArgs(trailingOnly = TRUE)

# Input
seurat_file  <- args[1]   # seurat_annotated_landscaper.RData
bindata_file <- args[2]   # BIN_DATA.tsv
group_file   <- args[3]   # group.tsv
allstates_file <- args[4] # Allstates.tsv
e_file       <- args[5]   # E.tsv
basin_file   <- args[6]   # Basin.tsv
mg_file      <- args[7]   # major_group.tsv
coord_file   <- args[8]   # Coordinate.tsv
igraph_file  <- args[9]   # igraph.RData
eb_file      <- args[10]  # EnergyBarrier.tsv
asmg_file    <- args[11]  # Allstates_major_group.tsv

# Output
output_file  <- args[12]

##----------------------------------------------------------------------
## Load
##----------------------------------------------------------------------
load(seurat_file)

BIN <- as.matrix(read.table(bindata_file, header = FALSE))

Allstates <- as.matrix(read.table(allstates_file, header = FALSE))
E <- as.numeric(read.table(e_file, header = FALSE)[, 1])
Basin <- as.vector(unlist(read.table(basin_file, header = FALSE)))
major_group_df <- read.table(mg_file, header = FALSE)
Coordinate <- as.matrix(read.table(coord_file, header = FALSE))
load(igraph_file)  # loads g
EnergyBarrier <- as.matrix(read.table(eb_file, header = FALSE))
Allstates_mg <- read.table(asmg_file, header = FALSE, stringsAsFactors = FALSE)

group <- read.table(group_file, header = FALSE)[, 1]

##----------------------------------------------------------------------
## Parse Allstates_major_group
##----------------------------------------------------------------------
K <- ncol(Allstates_mg) - 2
state_id_all <- Allstates_mg[, K + 1]
mg_label_all <- as.character(Allstates_mg[, K + 2])

##----------------------------------------------------------------------
## Compute state_idx: cell -> row index in Allstates
##----------------------------------------------------------------------
P <- nrow(Allstates)
key_fun <- function(v) paste0(v, collapse = "|")

keys_all  <- apply(Allstates, 1, key_fun)
keys_cell <- apply(BIN, 1, key_fun)
idx_map   <- setNames(seq_len(P), keys_all)
state_idx <- unname(idx_map[keys_cell])

if (any(is.na(state_idx))) {
    warning("Some cells have binary patterns not found in Allstates.")
}

##----------------------------------------------------------------------
## Align Allstates_major_group with Allstates
##----------------------------------------------------------------------
keys_mg  <- apply(as.matrix(Allstates_mg[, 1:K]), 1, key_fun)
mg_order <- match(keys_all, keys_mg)
if (any(is.na(mg_order))) {
    warning("Allstates_major_group rows do not fully match Allstates.")
}
state_id_aligned <- state_id_all[mg_order]
mg_label_aligned <- mg_label_all[mg_order]

##----------------------------------------------------------------------
## Cell-level meta.data
##----------------------------------------------------------------------
seurat.integrated$state_idx   <- state_idx
seurat.integrated$energy      <- E[state_idx]
seurat.integrated$is_basin    <- state_idx %in% Basin
seurat.integrated$state_id    <- state_id_aligned[state_idx]
seurat.integrated$major_group <- mg_label_aligned[state_idx]
seurat.integrated$landscaper_group <- group

##----------------------------------------------------------------------
## State-level data in @misc$landscaper
##----------------------------------------------------------------------
landscaper_data <- list(
    Allstates             = Allstates,
    E                     = E,
    Basin                 = Basin,
    major_group           = major_group_df,
    Coordinate            = Coordinate,
    igraph                = g,
    EnergyBarrier         = EnergyBarrier,
    Allstates_major_group = Allstates_mg,
    state_id              = state_id_aligned,
    mg_label              = mg_label_aligned,
    BIN_DATA              = BIN
)

seurat.integrated@misc$landscaper <- landscaper_data

##----------------------------------------------------------------------
## Save
##----------------------------------------------------------------------
save(seurat.integrated, file = output_file)
