source("src/Functions.R")

# Parameter
pc <- commandArgs(trailingOnly=TRUE)[1]
outfile1 <- commandArgs(trailingOnly=TRUE)[2]
outfile2 <- commandArgs(trailingOnly=TRUE)[3]

# Loading
all_states_integrated <- unlist(read.delim(paste0('plot/integrated/', pc, '/Landscaper/Allstates.tsv'), header=FALSE, sep="|"))
bin_data_integrated <- unlist(read.delim(paste0('output/integrated/', pc, '/binpca/BIN_DATA.tsv'), header=FALSE, sep="|"))
load('output/integrated/seurat_annotated_landscaper.RData')

# Sort
names(all_states_integrated) <- NULL
names(bin_data_integrated) <- NULL

# rm \t
bin_data_integrated <- gsub("\t", " ", bin_data_integrated)

# 各データごとの状態No
target_integrated <- match(bin_data_integrated, all_states_integrated)

# Assign Labels
names(target_integrated) <- colnames(seurat.integrated)
seurat.integrated$states <- target_integrated

# Plot
g <- DimPlot(seurat.integrated, reduction = "umap", group.by="states", label=FALSE, pt.size=1, label.size=6) + NoLegend() +
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
ggsave(file=outfile1, g, dpi=200, width=6, height=6)

# Plot
g <- DimPlot(seurat.integrated, reduction = "umap", group.by="states", split.by="sample", label=FALSE, ncol=3, pt.size=1, label.size=6) + NoLegend() +
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
ggsave(file=outfile2, g, dpi=200, width=20, height=13)
