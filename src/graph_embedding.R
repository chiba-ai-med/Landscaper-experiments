source("src/Functions.R")

# Parameter
args <- commandArgs(trailingOnly=TRUE)
infile1 <- args[1]  # P_metropolis.tsv
infile2 <- args[2]  # P_glauber.tsv
outfile1 <- args[3]
outfile2 <- args[4]

# Loading
P_m <- as.matrix(read.table(infile1, header=FALSE))
P_g <- as.matrix(read.table(infile2, header=FALSE))

# SVD
result_m <- svd(P_m)
result_g <- svd(P_g)

# Save
save(result_m, file=outfile1)
save(result_g, file=outfile2)
