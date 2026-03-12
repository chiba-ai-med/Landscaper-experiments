source("src/Functions.R")

# Parameter
pc <- commandArgs(trailingOnly=TRUE)[1]
outfile1 <- commandArgs(trailingOnly=TRUE)[2]
outfile2 <- commandArgs(trailingOnly=TRUE)[3]

# Loading
h_cont <- read.table(paste0("plot/cont/", pc, "/Landscaper/h.tsv"), header=FALSE)
h_DAPT <- read.table(paste0("plot/DAPT/", pc, "/Landscaper/h.tsv"), header=FALSE)
h_cont_cov <- read.table(paste0("plot/cont_cov/", pc, "/Landscaper/h.tsv"), header=FALSE)
h_DAPT_cov <- read.table(paste0("plot/DAPT_cov/", pc, "/Landscaper/h.tsv"), header=FALSE)

# Name (w/o cov)
rownames_h <- paste0("Var", seq(nrow(h_cont)))

# Re-ordering (w/o cov)
h_diff <- h_DAPT - h_cont
order_h_diff <- order(unlist(h_diff), decreasing=TRUE)
h_cont <- h_cont[order_h_diff, ]
h_DAPT <- h_DAPT[order_h_diff, ]
h_diff <- h_diff[order_h_diff, ]
rownames_h <- rownames_h[order_h_diff]

# Name (w/ cov)
rownames_h_cov <- paste0("Var", seq(nrow(h_cont_cov)))

# Re-ordering (w/ cov)
h_diff_cov <- h_DAPT_cov - h_cont_cov
order_h_diff_cov <- order(unlist(h_diff_cov), decreasing=TRUE)
h_cont_cov <- h_cont_cov[order_h_diff_cov, ]
h_DAPT_cov <- h_DAPT_cov[order_h_diff_cov, ]
h_diff_cov <- h_diff_cov[order_h_diff_cov, ]
rownames_h_cov <- rownames_h_cov[order_h_diff_cov]

# Data for ggplot2 (w/o cov)
gdata <- data.frame(value=c(h_cont, h_DAPT),
        group=c(rep("cont", length=length(h_cont)),
            rep("DAPT", length=length(h_DAPT))),
        variable=factor(rownames_h, level=rownames_h))

# Data for ggplot2 (w/ cov)
gdata_cov <- data.frame(value=c(h_cont_cov, h_DAPT_cov),
        group=c(rep("cont", length=length(h_cont_cov)),
            rep("DAPT", length=length(h_DAPT_cov))),
        variable=factor(rownames_h_cov, level=rownames_h_cov))

# Setting
min_val <- min(gdata$value, gdata_cov$value)
max_val <- max(gdata$value, gdata_cov$value)

# Plot (w/o cov)
g1 <- ggplot() + geom_line(aes(x=variable, y=value, colour=group, group=group), linewidth=1, data = gdata, stat="identity") + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + ylim(min_val, max_val)
ggsave(outfile1, plot=g1, width=7, height=4.5)

# Plot (w cov)
g2 <- ggplot() + 
    geom_line(aes(x=variable, y=value, colour=group, group=group), linewidth=1, data = gdata_cov, stat="identity") + 
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
    ylim(min_val, max_val)
ggsave(outfile2, plot=g2, width=7, height=4.5)