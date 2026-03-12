source("src/Functions.R")

# Parameter
pc <- commandArgs(trailingOnly=TRUE)[1]
outfile1 <- commandArgs(trailingOnly=TRUE)[2]
outfile2 <- commandArgs(trailingOnly=TRUE)[3]

# Loading
all_states_cont_cov <- unlist(read.delim(paste0('plot/cont_cov/', pc, '/Landscaper/Allstates.tsv'), header=FALSE, sep="|"))
bin_data_cont_cov <- unlist(read.delim(paste0('output/cont_cov/', pc, '/binpca/BIN_DATA.tsv'), header=FALSE, sep="|"))
load('output/cont_cov/seurat_annotated_landscaper.RData')

# Sort
names(all_states_cont_cov) <- NULL
names(bin_data_cont_cov) <- NULL

# 各データごとの状態No
target_cont_cov <- match(bin_data_cont_cov, all_states_cont_cov)

# Assign Labels
names(target_cont_cov) <- colnames(seurat.integrated)
seurat.integrated$states <- target_cont_cov

# Plot
g <- DimPlot(seurat.integrated, reduction = "umap", group.by="states", label=FALSE, pt.size=2, label.size=6) + NoLegend() +
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
png(file=outfile1, width=600, height=600)
print(g)
dev.off()

# Plot
g <- DimPlot(seurat.integrated, reduction = "umap", group.by="states", split.by="sample", label=FALSE, pt.size=2, label.size=6) + NoLegend() +
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
png(file=outfile2, width=2000, height=650)
print(g)
dev.off()

