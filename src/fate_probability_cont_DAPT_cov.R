source("src/Functions.R")

# Parameter
args <- commandArgs(trailingOnly=TRUE)
infile1 <- args[1]  # cont_cov seurat_landscaper.RData (for Basin)
infile2 <- args[2]  # DAPT_cov seurat_landscaper.RData (for Basin)
infile3 <- args[3]  # cont_cov P_metropolis.tsv
infile4 <- args[4]  # DAPT_cov P_metropolis.tsv
infile5 <- args[5]  # cont_cov P_glauber.tsv
infile6 <- args[6]  # DAPT_cov P_glauber.tsv
outfile1 <- args[7]
outfile2 <- args[8]
outfile3 <- args[9]
outfile4 <- args[10]

# Loading
cont_env <- new.env()
load(infile1, envir = cont_env)
Basin_cont <- cont_env$seurat.integrated@misc$landscaper$Basin

dapt_env <- new.env()
load(infile2, envir = dapt_env)
Basin_DAPT <- dapt_env$seurat.integrated@misc$landscaper$Basin
rm(cont_env, dapt_env)
Basin <- union(Basin_cont, Basin_DAPT)

# cont
P_m <- as.matrix(read.table(infile3, header=FALSE))
P_g <- as.matrix(read.table(infile5, header=FALSE))

## Metropolis
res_m <- absorption_probabilities(P_m, absorbing = Basin, reg = 1e-4)
H_m <- fate_entropy(res_m$F, base = 2)
argmax_m <- fate_argmax(res_m$F)

## Glauber
res_g <- absorption_probabilities(P_g, absorbing = Basin, reg = 1e-4)
H_g <- fate_entropy(res_g$F, base = 2)
argmax_g <- fate_argmax(res_g$F)

## Save
save(res_m, H_m, argmax_m, file=outfile1)
save(res_g, H_g, argmax_g, file=outfile2)


# DAPT
P_m <- as.matrix(read.table(infile4, header=FALSE))
P_g <- as.matrix(read.table(infile5, header=FALSE))

## Metropolis
res_m <- absorption_probabilities(P_m, absorbing = Basin, reg = 1e-4)
H_m <- fate_entropy(res_m$F, base = 2)
argmax_m <- fate_argmax(res_m$F)

## Glauber
res_g <- absorption_probabilities(P_g, absorbing = Basin, reg = 1e-4)
H_g <- fate_entropy(res_g$F, base = 2)
argmax_g <- fate_argmax(res_g$F)

## Save
save(res_m, H_m, argmax_m, file=outfile3)
save(res_g, H_g, argmax_g, file=outfile4)
