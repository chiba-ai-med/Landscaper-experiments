source("src/Functions.R")

# Parameter
args <- commandArgs(trailingOnly=TRUE)
infile1 <- args[1]  # seurat_landscaper.RData (for Basin)
infile2 <- args[2]  # P_metropolis.tsv
infile3 <- args[3]  # P_glauber.tsv
outfile1 <- args[4]
outfile2 <- args[5]

# Loading
load(infile1)
Basin <- seurat.integrated@misc$landscaper$Basin
P_m <- as.matrix(read.table(infile2, header=FALSE))
P_g <- as.matrix(read.table(infile3, header=FALSE))

# Metropolis
res_m <- absorption_probabilities(P_m, absorbing = Basin, reg = 1e-4)
H_m <- fate_entropy(res_m$F, base = 2)
argmax_m <- fate_argmax(res_m$F)

# Glauber
res_g <- absorption_probabilities(P_g, absorbing = Basin, reg = 1e-4)
H_g <- fate_entropy(res_g$F, base = 2)
argmax_g <- fate_argmax(res_g$F)

# Save
save(res_m, H_m, argmax_m, file=outfile1)
save(res_g, H_g, argmax_g, file=outfile2)
