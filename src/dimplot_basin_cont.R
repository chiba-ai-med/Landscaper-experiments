source("src/Functions.R")

# Parameter
pc <- commandArgs(trailingOnly=TRUE)[1]
outfile1 <- commandArgs(trailingOnly=TRUE)[2]
outfile2 <- commandArgs(trailingOnly=TRUE)[3]

# Loading
all_states_cont <- unlist(read.delim(paste0('plot/cont/', pc, '/Landscaper/Allstates.tsv'), header=FALSE, sep="|"))
bin_data_cont <- unlist(read.delim(paste0('output/cont/', pc, '/binpca/BIN_DATA.tsv'), header=FALSE, sep="|"))
basin_cont <- unlist(read.table(paste0('plot/cont/', pc, '/Landscaper/Basin.tsv'), header=FALSE))
load('output/cont/seurat_annotated_landscaper.RData')

# Sort
names(all_states_cont) <- NULL
names(bin_data_cont) <- NULL

# 各データごとの状態No
target_cont <- match(bin_data_cont, all_states_cont)

basin_cont_sorted <- rep(0, length=length(bin_data_cont))

target_cont2 <- unlist(sapply(basin_cont, function(x){
     which(target_cont == x)
}))

basin_cont_sorted[target_cont2] <- 1

# Assign Labels
seurat.integrated$basin <- basin_cont_sorted

# Plot
g <- DimPlot(seurat.integrated, reduction = "umap", group.by="basin", label=FALSE, pt.size=2, label.size=6, cols=c(8,3)) + NoLegend() +
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
g <- DimPlot(seurat.integrated, reduction = "umap", group.by="basin", split.by="sample", label=FALSE, pt.size=2, label.size=6, cols=c(8,3)) + NoLegend() +
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
png(file=outfile2, width=2400, height=650)
print(g)
dev.off()
