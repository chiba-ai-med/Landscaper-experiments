source("src/Functions.R")

# Parameter
pc <- commandArgs(trailingOnly=TRUE)[1]
outfile <- commandArgs(trailingOnly=TRUE)[2]

# Loading
all_states_integrated <- unlist(read.delim(paste0('plot/oulhen/integrated/', pc, '/Landscaper/Allstates.tsv'), header=FALSE, sep="|"))
bin_data_integrated <- unlist(read.delim(paste0('output/oulhen/integrated/', pc, '/binpca/BIN_DATA.tsv'), header=FALSE, sep="|"))
basin_integrated <- unlist(read.table(paste0('plot/oulhen/integrated/', pc, '/Landscaper/Basin.tsv'), header=FALSE))
load('output/oulhen/integrated/seurat_annotated_landscaper.RData')

# Sort
names(all_states_integrated) <- NULL
names(bin_data_integrated) <- NULL

bin_data_integrated <- gsub("\t", " ", bin_data_integrated)

# 各データごとの状態No
target_integrated <- match(bin_data_integrated, all_states_integrated)

basin_integrated_sorted <- rep(0, length=length(bin_data_integrated))

target_integrated2 <- unlist(sapply(basin_integrated, function(x){
     which(target_integrated == x)
}))

basin_integrated_sorted[target_integrated2] <- 1

# Assign Labels
seurat.integrated$basin <- basin_integrated_sorted

# Plot
g <- DimPlot(seurat.integrated, reduction = "umap", group.by="basin", label=FALSE, pt.size=3, label.size=6, cols=c(8,3)) + NoLegend() +
theme(axis.line = element_blank(),
           axis.text.x = element_blank(),
           axis.text.y = element_blank(),
           axis.ticks = element_blank(),
           axis.title.x = element_blank(),
           axis.title.y = element_blank(),
           panel.background = element_blank(),
           panel.border = element_blank(),
           panel.grid.major = element_blank(),
           panel.grid.minor = element_blank(),
           plot.background = element_blank())
ggsave(file=outfile, g, dpi=200, width=6, height=6)
