library("akima")
library("igraph")

args <- commandArgs(trailingOnly = TRUE)
pc <- args[1]
outfile1 <- args[2]
outfile2 <- args[3]
outfile3 <- args[4]
outfile4 <- args[5]

# Load
E_cont <- unlist(read.table(paste0("plot/cont/", pc, "/Landscaper/E.tsv"), header=FALSE))
E_DAPT <- unlist(read.table(paste0("plot/DAPT/", pc, "/Landscaper/E.tsv"), header=FALSE))
E_cont_cov <- unlist(read.table(paste0("plot/cont_cov/", pc, "/Landscaper/E.tsv"), header=FALSE))
E_DAPT_cov <- unlist(read.table(paste0("plot/DAPT_cov/", pc, "/Landscaper/E.tsv"), header=FALSE))

G_sub_cont <- unlist(read.table(paste0("plot/cont/", pc, "/Landscaper/SubGraph.tsv"), header=FALSE))
G_sub_DAPT <- unlist(read.table(paste0("plot/DAPT/", pc, "/Landscaper/SubGraph.tsv"), header=FALSE))
G_sub_cont_cov <- unlist(read.table(paste0("plot/cont_cov/", pc, "/Landscaper/SubGraph.tsv"), header=FALSE))
G_sub_DAPT_cov <- unlist(read.table(paste0("plot/DAPT_cov/", pc, "/Landscaper/SubGraph.tsv"), header=FALSE))

Basin_cont <- unlist(read.table(paste0("plot/cont/", pc, "/Landscaper/Basin.tsv"), header=FALSE))
Basin_DAPT <- unlist(read.table(paste0("plot/DAPT/", pc, "/Landscaper/Basin.tsv"), header=FALSE))
Basin_cont_cov <- unlist(read.table(paste0("plot/cont_cov/", pc, "/Landscaper/Basin.tsv"), header=FALSE))
Basin_DAPT_cov <- unlist(read.table(paste0("plot/DAPT_cov/", pc, "/Landscaper/Basin.tsv"), header=FALSE))

Coordinate <- as.matrix(read.table(paste0("plot/integrated/", pc, "/Landscaper/Coordinate.tsv"), header=FALSE))
load(paste0("plot/integrated/", pc, "/Landscaper/igraph.RData"))

# Basin強調用のサイズ・色ベクトルを生成するヘルパー
n_states <- length(V(g))
make_basin_style <- function(basin_idx) {
  vcol <- rep(rgb(0,0,0), n_states)
  vcol[basin_idx] <- rgb(1,0,0)
  vsize <- rep(3, n_states)
  vsize[basin_idx] <- 6
  lcex <- rep(3.0, n_states)
  lcex[basin_idx] <- 6.0
  list(color = vcol, size = vsize, lcex = lcex)
}

style_cont <- make_basin_style(Basin_cont)
style_DAPT <- make_basin_style(Basin_DAPT)
style_cont_cov <- make_basin_style(Basin_cont_cov)
style_DAPT_cov <- make_basin_style(Basin_DAPT_cov)

# Preprocess
tmp_cont <- interp(Coordinate[, 1], Coordinate[, 2], E_cont, nx = 200, ny = 200)
tmp_DAPT <- interp(Coordinate[, 1], Coordinate[, 2], E_DAPT, nx = 200, ny = 200)
tmp_cont_cov <- interp(Coordinate[, 1], Coordinate[, 2], E_cont_cov, nx = 200, ny = 200)
tmp_DAPT_cov <- interp(Coordinate[, 1], Coordinate[, 2], E_DAPT_cov, nx = 200, ny = 200)

zlim_cont <- tmp_cont$z |> as.vector() |> na.omit() |> range()
zlim_DAPT <- tmp_DAPT$z |> as.vector() |> na.omit() |> range()
zlim_cont_cov <- tmp_cont_cov$z |> as.vector() |> na.omit() |> range()
zlim_DAPT_cov <- tmp_DAPT_cov$z |> as.vector() |> na.omit() |> range()
zlim <- range(c(zlim_cont, zlim_DAPT, zlim_cont_cov, zlim_DAPT_cov))

# Plot
png(file=outfile1, width=1500, height=1500)
filled.contour(tmp_cont, color.palette = colorRampPalette(topo.colors(11, alpha = 1)),
    xlim = c(-1, 1), ylim = c(-1, 1), zlim = zlim, asp = 1, nlevels = 24,
    key.title = title(main = paste0("Energy\n[", round(zlim[1], 2),
                                    ", ", round(zlim[2], 2), "]"), cex.main = 1.5),
    axes = FALSE,
    plot.axes = {
    contour(tmp_cont$x, tmp_cont$y, tmp_cont$z, add = TRUE,
        levels = seq(-10, 10, by=1), col = "gray75")
    plot(g, layout=Coordinate,
        vertex.frame.color = NA,
        vertex.label.color = style_cont$color,
        vertex.label.cex = style_cont$lcex,
        vertex.size = 0.2 * style_cont$size,
        vertex.color = factor(G_sub_cont),
        edge.arrow.size = 0.2, edge.width = 2,
        add = TRUE, axes=FALSE)
    })
dev.off()

png(file=outfile2, width=1500, height=1500)
filled.contour(tmp_DAPT, color.palette = colorRampPalette(topo.colors(11, alpha = 1)),
    xlim = c(-1, 1), ylim = c(-1, 1), zlim = zlim, asp = 1, nlevels = 24,
    key.title = title(main = paste0("Energy\n[", round(zlim[1], 2),
                                    ", ", round(zlim[2], 2), "]"), cex.main = 1.5),
    axes = FALSE,
    plot.axes = {
    contour(tmp_DAPT$x, tmp_DAPT$y, tmp_DAPT$z, add = TRUE,
        levels = seq(-10, 10, by=1), col = "gray75")
    plot(g, layout=Coordinate,
        vertex.frame.color = NA,
        vertex.label.color = style_DAPT$color,
        vertex.label.cex = style_DAPT$lcex,
        vertex.size = 0.2 * style_DAPT$size,
        vertex.color = factor(G_sub_DAPT),
        edge.arrow.size = 0.2, edge.width = 2,
        add = TRUE, axes=FALSE)
    })
dev.off()

png(file=outfile3, width=1500, height=1500)
filled.contour(tmp_cont_cov, color.palette = colorRampPalette(topo.colors(11, alpha = 1)),
    xlim = c(-1, 1), ylim = c(-1, 1), zlim = zlim, asp = 1, nlevels = 24,
    key.title = title(main = paste0("Energy\n[", round(zlim[1], 2),
                                    ", ", round(zlim[2], 2), "]"), cex.main = 1.5),
    axes = FALSE,
    plot.axes = {
    contour(tmp_cont_cov$x, tmp_cont_cov$y, tmp_cont_cov$z, add = TRUE,
        levels = seq(-10, 10, by=1), col = "gray75")
    plot(g, layout=Coordinate,
        vertex.frame.color = NA,
        vertex.label.color = style_cont_cov$color,
        vertex.label.cex = style_cont_cov$lcex,
        vertex.size = 0.2 * style_cont_cov$size,
        vertex.color = factor(G_sub_cont_cov),
        edge.arrow.size = 0.2, edge.width = 2,
        add = TRUE, axes=FALSE)
    })
dev.off()

png(file=outfile4, width=1500, height=1500)
filled.contour(tmp_DAPT_cov, color.palette = colorRampPalette(topo.colors(11, alpha = 1)),
    xlim = c(-1, 1), ylim = c(-1, 1), zlim = zlim, asp = 1, nlevels = 24,
    key.title = title(main = paste0("Energy\n[", round(zlim[1], 2),
                                    ", ", round(zlim[2], 2), "]"), cex.main = 1.5),
    axes = FALSE,
    plot.axes = {
    contour(tmp_DAPT_cov$x, tmp_DAPT_cov$y, tmp_DAPT_cov$z, add = TRUE,
        levels = seq(-10, 10, by=1), col = "gray75")
    plot(g, layout=Coordinate,
        vertex.frame.color = NA,
        vertex.label.color = style_DAPT_cov$color,
        vertex.label.cex = style_DAPT_cov$lcex,
        vertex.size = 0.2 * style_DAPT_cov$size,
        vertex.color = factor(G_sub_DAPT_cov),
        edge.arrow.size = 0.2, edge.width = 2,
        add = TRUE, axes=FALSE)
    })
dev.off()
