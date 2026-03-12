args <- commandArgs(trailingOnly=TRUE)
outfile <- args[length(args)]
infiles <- args[-length(args)]

vals <- unlist(lapply(infiles, function(f) unlist(read.table(f, header=FALSE))))
write.table(
  data.frame(min=min(vals), max=max(vals)),
  outfile, row.names=FALSE, col.names=TRUE, quote=FALSE, sep="\t"
)
