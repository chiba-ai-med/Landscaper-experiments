source("src/Functions.R")

# Parameter
pc <- commandArgs(trailingOnly=TRUE)[1]
outfile <- commandArgs(trailingOnly=TRUE)[2]

# Loading
all_states_DAPT <- unlist(read.delim(paste0('plot/oulhen/DAPT/', pc, '/Landscaper/Allstates.tsv'), header=FALSE, sep="|"))
bin_data_DAPT <- unlist(read.delim(paste0('output/oulhen/DAPT/', pc, '/binpca/BIN_DATA.tsv'), header=FALSE, sep="|"))
load('output/oulhen/DAPT/seurat_annotated_landscaper.RData')

# Sort
names(all_states_DAPT) <- NULL
names(bin_data_DAPT) <- NULL

bin_data_DAPT <- gsub("\t", " ", bin_data_DAPT)

# 各データごとの状態No
target_DAPT <- match(bin_data_DAPT, all_states_DAPT)

# Assign Labels
names(target_DAPT) <- colnames(seurat.integrated)
seurat.integrated$states <- target_DAPT

# Plot
g <- DimPlot(seurat.integrated, reduction = "umap", group.by="states", label=FALSE, pt.size=3, label.size=6) + NoLegend() +
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
png(file=outfile, width=600, height=600)
print(g)
dev.off()
