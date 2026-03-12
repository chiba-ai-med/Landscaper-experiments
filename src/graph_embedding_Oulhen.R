source("src/Functions.R")

# Parameter
infile1 <- commandArgs(trailingOnly=TRUE)[1]
infile2 <- commandArgs(trailingOnly=TRUE)[2]
infile3 <- commandArgs(trailingOnly=TRUE)[3]
outfile1 <- commandArgs(trailingOnly=TRUE)[4]
outfile2 <- commandArgs(trailingOnly=TRUE)[5]
# infile1 <- 'plot/hpbase/integrated/Landscaper/Allstates.tsv'
# infile2 <- 'plot/hpbase/integrated/P_metropolis.tsv'
# infile3 <- 'plot/hpbase/integrated/P_glauber.tsv'

# Loading
Allstates <- as.matrix(read.table(infile1, header=FALSE))
P_m <- as.matrix(read.table(infile2, header=FALSE))
P_g <- as.matrix(read.table(infile3, header=FALSE))

# SVD
result_m <- svd(P_m)
result_g <- svd(P_g)

# Save
save(result_m, file=outfile1)
save(result_g, file=outfile2)
