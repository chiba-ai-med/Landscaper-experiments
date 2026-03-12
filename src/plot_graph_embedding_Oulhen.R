source("src/Functions.R")

# Parameter
infile1 <- commandArgs(trailingOnly=TRUE)[1]
infile2 <- commandArgs(trailingOnly=TRUE)[2]
infile3 <- commandArgs(trailingOnly=TRUE)[3]
infile4 <- commandArgs(trailingOnly=TRUE)[4]
infile5 <- commandArgs(trailingOnly=TRUE)[5]
outfile1 <- commandArgs(trailingOnly=TRUE)[6]
outfile2 <- commandArgs(trailingOnly=TRUE)[7]
# infile1 <- 'plot/hpbase/integrated/Landscaper/E.tsv'
# infile2 <- 'plot/hpbase/integrated/Landscaper/Basin.tsv'
# infile3 <- 'plot/hpbase/integrated/Landscaper/igraph.RData'
# infile4 <- 'plot/hpbase/integrated/P_metropolis_emded.RData'
# infile5 <- 'plot/hpbase/integrated/P_glauber_emded.RData'

# Load
E <- unlist(read.table(infile1, header=FALSE))
Basin <- unlist(read.table(infile2, header=FALSE))
load(infile3)
load(infile4)
load(infile5)

# ベイシンの設定
basin_color <- rep(rgb(0,0,0), length=length(V(g)))
basin_color[Basin] <- rgb(1,0,0)
basin_size <- rep(0.75, length=length(V(g)))
basin_size[Basin] <- 1.5

# Plot (metropolis)
## 座標（2次元）
nc <- ncol(result_m$u)
coords <- result_m$u[, 1:2, drop = FALSE]

# 1. Kamada–Kawai で初期位置を少し調整
lay_init <- layout_with_kk(g, coords = coords)

# 2. FR layout で “拡散”
set.seed(1234)
lay <- layout_with_fr(g,
                      coords = lay_init,     # 初期座標を指定
                      niter = 300,
                      grid = "nogrid",
                      weights = E(g)$weight) # 重みを考慮してもOK

# 3. 正規化（±1 範囲にスケール）
coords <- layout.norm(lay, -1, 1, -1, 1)

## 補間（akima::interp 推奨）
## install.packages("akima") が必要な場合あり
xo <- seq(min(coords[,1]), max(coords[,1]), length.out = 200)
yo <- seq(min(coords[,2]), max(coords[,2]), length.out = 200)
tmp <- akima::interp(x = coords[,1], y = coords[,2], z = E,
                     xo = xo, yo = yo, duplicate = "mean", linear = TRUE)

## 色レンジ
zlim <- range(tmp$z, finite = TRUE)

png(file = outfile1, width = 1500, height = 1500)
filled.contour(x = tmp$x, y = tmp$y, z = tmp$z,
  color.palette = colorRampPalette(topo.colors(11)),
  xlim = range(tmp$x), ylim = range(tmp$y),
  zlim = zlim, asp = 1, nlevels = 24,
  key.title = title(main = "Energy"),
  axes = FALSE,
  plot.axes = {
    contour(tmp$x, tmp$y, tmp$z, add = TRUE,
            levels = pretty(zlim, n = 20), col = "gray75")
    # ここが重要：rescale=FALSE で既存座標系にそのまま描く
    plot(g,
         layout = coords,
         rescale = FALSE,
         xlim = range(tmp$x), ylim = range(tmp$y),
         vertex.frame.color = NA,
         vertex.label.color = basin_color,
         vertex.label.cex = basin_size,
         vertex.size = 0.1*basin_size,
         edge.arrow.size = 0.1,
         add = TRUE, axes = FALSE)
  })
dev.off()


# Plot (glauber)
## 座標（2次元）
nc <- ncol(result_g$u)
coords <- result_g$u[, 1:2, drop = FALSE]

# 1. Kamada–Kawai で初期位置を少し調整
lay_init <- layout_with_kk(g, coords = coords)

# 2. FR layout で “拡散”
set.seed(1234)
lay <- layout_with_fr(g,
                      coords = lay_init,     # 初期座標を指定
                      niter = 300,
                      grid = "nogrid",
                      weights = E(g)$weight) # 重みを考慮してもOK

# 3. 正規化（±1 範囲にスケール）
coords <- layout.norm(lay, -1, 1, -1, 1)

## 補間（akima::interp 推奨）
## install.packages("akima") が必要な場合あり
xo <- seq(min(coords[,1]), max(coords[,1]), length.out = 200)
yo <- seq(min(coords[,2]), max(coords[,2]), length.out = 200)
tmp <- akima::interp(x = coords[,1], y = coords[,2], z = E,
                     xo = xo, yo = yo, duplicate = "mean", linear = TRUE)

## 色レンジ
zlim <- range(tmp$z, finite = TRUE)

png(file = outfile2, width = 1500, height = 1500)
filled.contour(x = tmp$x, y = tmp$y, z = tmp$z,
  color.palette = colorRampPalette(topo.colors(11)),
  xlim = range(tmp$x), ylim = range(tmp$y),
  zlim = zlim, asp = 1, nlevels = 24,
  key.title = title(main = "Energy"),
  axes = FALSE,
  plot.axes = {
    contour(tmp$x, tmp$y, tmp$z, add = TRUE,
            levels = pretty(zlim, n = 20), col = "gray75")
    # ここが重要：rescale=FALSE で既存座標系にそのまま描く
    plot(g,
         layout = coords,
         rescale = FALSE,
         xlim = range(tmp$x), ylim = range(tmp$y),
         vertex.frame.color = NA,
         vertex.label.color = basin_color,
         vertex.label.cex = basin_size,
         vertex.size = 0.1*basin_size,
         edge.arrow.size = 0.1,
         add = TRUE, axes = FALSE)
  })
dev.off()

