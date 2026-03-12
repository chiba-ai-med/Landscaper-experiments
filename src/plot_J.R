source("src/Functions.R")

# Parameter
pc <- commandArgs(trailingOnly=TRUE)[1]
outfile1 <- commandArgs(trailingOnly=TRUE)[2]
outfile2 <- commandArgs(trailingOnly=TRUE)[3]

# Loading
J_cont <- read.table(paste0("plot/cont/", pc, "/Landscaper/J.tsv"), header=FALSE)
J_DAPT <- read.table(paste0("plot/DAPT/", pc, "/Landscaper/J.tsv"), header=FALSE)
J_cont_cov <- read.table(paste0("plot/cont_cov/", pc, "/Landscaper/J.tsv"), header=FALSE)
J_DAPT_cov <- read.table(paste0("plot/DAPT_cov/", pc, "/Landscaper/J.tsv"), header=FALSE)

# Name
upt_J <- upper.tri(J_cont)
varnames_J <- c()
for(i in seq(nrow(upt_J))){
	for(j in seq(ncol(upt_J))){
		if(upt_J[i,j]){
			varnames_J <- c(varnames_J, paste0("Var", i, "-Var", j))
		}
	}
}
varnames_J_cov <- varnames_J

# Matrix => Vector
J_cont <- J_cont[upt_J]
J_DAPT <- J_DAPT[upt_J]
J_cont_cov <- J_cont_cov[upt_J]
J_DAPT_cov <- J_DAPT_cov[upt_J]

# Re-ordering (w/o cov)
J_diff <- J_DAPT - J_cont
order_J_diff <- order(unlist(J_diff), decreasing=TRUE)
J_cont <- J_cont[order_J_diff]
J_DAPT <- J_DAPT[order_J_diff]
J_diff <- J_diff[order_J_diff]
varnames_J <- varnames_J[order_J_diff]

# Re-ordering (w/ cov)
J_diff_cov <- J_DAPT_cov - J_cont_cov
order_J_diff_cov <- order(unlist(J_diff_cov), decreasing=TRUE)
J_cont_cov <- J_cont_cov[order_J_diff_cov]
J_DAPT_cov <- J_DAPT_cov[order_J_diff_cov]
J_diff_cov <- J_diff_cov[order_J_diff_cov]
varnames_J_cov <- varnames_J_cov[order_J_diff_cov]

# Data for ggplot2 (w/o cov)
gdata <- data.frame(value=c(J_cont, J_DAPT),
        group=c(rep("cont", length=length(J_cont)),
            rep("DAPT", length=length(J_DAPT))),
        variable=factor(varnames_J, level=varnames_J))

# Data for ggplot2 (w/ cov)
gdata_cov <- data.frame(value=c(J_cont_cov, J_DAPT_cov),
        group=c(rep("cont", length=length(J_cont_cov)),
            rep("DAPT", length=length(J_DAPT_cov))),
        variable=factor(varnames_J_cov, level=varnames_J_cov))

# Setting
min_val <- min(gdata$value, gdata_cov$value)
max_val <- max(gdata$value, gdata_cov$value)

# Plot (w/o cov)
g1 <- ggplot() + geom_line(aes(x=variable, y=value, colour=group, group=group), linewidth=1, data = gdata, stat="identity") + theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(outfile1, plot=g1, width=12, height=4.5)

# Plot (w cov)
g2 <- ggplot() + geom_line(aes(x=variable, y=value, colour=group, group=group), linewidth=1, data = gdata_cov, stat="identity") + theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(outfile2, plot=g2, width=12, height=4.5)
