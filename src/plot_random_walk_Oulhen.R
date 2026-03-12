source("src/Functions.R")

# Parameter
args <- commandArgs(trailingOnly=TRUE)
infile1  <- args[1]  # major_group.tsv
infile2  <- args[2]  # random_walk list output (.RData)
outfile1 <- args[3]  # Aboral
outfile2 <- args[4]  # Oral
outfile3 <- args[5]  # Ciliary
outfile4 <- args[6]  # Neurons

# Load labels
major_group <- read.table(infile1, header=FALSE, sep="", quote="",
                          stringsAsFactors=FALSE, check.names=FALSE)
labels_base <- trimws(sub("\\..*$", "", as.character(major_group[[ncol(major_group)]])))

# Load P_list (RData)
load(infile2)
if (exists("P_g_list")) {
  P_list <- P_g_list
} else if (exists("P_m_list")) {
  P_list <- P_m_list
} else {
  stop("P_g_list または P_m_list が見つかりません")
}

# starts and destinations
start_levels <- c("Aboral_ectoderm","Oral_ectoderm","Ciliary_band","Neurons")
dest_levels  <- c("Aboral_ectoderm","Oral_ectoderm","Ciliary_band","Neurons","Other")
focus_types  <- c("Aboral_ectoderm","Oral_ectoderm","Ciliary_band","Neurons")
celltype2 <- ifelse(labels_base %in% focus_types, labels_base, "Other")

# steps (0..n)
steps <- 0:(length(P_list) - 1)

# colors
pal <- c(
  "Aboral_ectoderm" = "#008080",
  "Oral_ectoderm"   = "#FFFF00",
  "Ciliary_band"    = "#00FF26",
  "Neurons"         = "#FF00FF",
  "Other"           = "grey70"
)

plot_one_start <- function(start_label, outfile) {
  idx_start <- which(labels_base == start_label)

  # If start label not found, output an "empty" plot (all zeros)
  if (length(idx_start) == 0) {
    df0 <- expand.grid(
      step = steps,
      dest = dest_levels,
      KEEP.OUT.ATTRS = FALSE,
      stringsAsFactors = FALSE
    )
    df0$prob <- 0
    df0$dest <- factor(df0$dest, levels = dest_levels)

    p <- ggplot(df0, aes(x = step, y = prob, color = dest)) +
      geom_line(linewidth = 2.0) +
      geom_point(size = 2.5) +
      scale_color_manual(values = pal, guide = "none") +
      scale_x_continuous(breaks = steps) +
      coord_cartesian(ylim = c(0, 1)) +
      labs(x = "Step", y = "Probability") +
      theme_classic(base_size = 12) +
      theme(text = element_text(size = 36))

    ggsave(outfile, p, width = 8, height = 6, dpi = 300)
    return(invisible(NULL))
  }

  df_list <- vector("list", length(P_list))
  for (k in seq_along(P_list)) {
    Pk <- P_list[[k]]
    P_dest_rows <- Pk[idx_start, , drop=FALSE]
    p_dest <- colMeans(P_dest_rows)

    tmp <- tapply(p_dest, celltype2, sum)
    v <- setNames(rep(0, length(dest_levels)), dest_levels)
    hit <- intersect(names(tmp), dest_levels)
    if (length(hit) > 0) v[hit] <- as.numeric(tmp[hit])

    df_list[[k]] <- tibble(
      step = steps[k],
      Aboral_ectoderm = v["Aboral_ectoderm"],
      Oral_ectoderm   = v["Oral_ectoderm"],
      Ciliary_band    = v["Ciliary_band"],
      Neurons         = v["Neurons"],
      Other           = v["Other"]
    )
  }

  df <- bind_rows(df_list) %>%
    pivot_longer(-step, names_to="dest", values_to="prob") %>%
    mutate(dest = factor(dest, levels = dest_levels))

  p <- ggplot(df, aes(x = step, y = prob, color = dest)) +
    geom_line(linewidth = 2.0) +
    geom_point(size = 2.5) +
    scale_color_manual(values = pal, guide = "none") +
    scale_x_continuous(breaks = steps) +
    coord_cartesian(ylim = c(0, 1)) +
    labs(x = "Step", y = "Probability") +
    theme_classic(base_size = 12) +
    theme(text = element_text(size = 36))

  ggsave(outfile, p, width = 8, height = 6, dpi = 300)
}

# run (5 starts -> 5 files)
plot_one_start("Aboral_ectoderm", outfile1)
plot_one_start("Oral_ectoderm",   outfile2)
plot_one_start("Ciliary_band",    outfile3)
plot_one_start("Neurons",         outfile4)