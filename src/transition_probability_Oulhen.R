source("src/Functions.R")

# Parameter
infile1 <- commandArgs(trailingOnly=TRUE)[1]
infile2 <- commandArgs(trailingOnly=TRUE)[2]
outfile1 <- commandArgs(trailingOnly=TRUE)[3]
outfile2 <- commandArgs(trailingOnly=TRUE)[4]
# infile1 <- 'plot/hpbase/integrated/Landscaper/E.tsv'
# infile2 <- 'plot/hpbase/integrated/Landscaper/Allstates.tsv'

##------------------------------------------------------------------------------
## ロード
##------------------------------------------------------------------------------

E         <- as.numeric(read.table(infile1, header = FALSE)[, 1])
Allstates <- as.matrix(read.table(infile2, header = FALSE))

P    <- nrow(Allstates)
Sdim <- ncol(Allstates)

##------------------------------------------------------------------------------
## 転移確率行列の構築（Metropolis / Glauber） with 時系列制約
##------------------------------------------------------------------------------

## 許容幅 tol:
## time_state[j] < time_state[i] - tol を禁止
## 36,48,72,96 なら tol=0 で「厳密に過去は禁止、同時刻は許可」
tol_h <- 0

P_metropolis <- transition_matrix_from_E(
  E          = E,
  S          = Allstates,
  beta       = 1,
  kernel     = "metropolis",
  tol        = tol_h
)

P_glauber <- transition_matrix_from_E(
  E          = E,
  S          = Allstates,
  beta       = 1,
  kernel     = "glauber",
  tol        = tol_h
)

##------------------------------------------------------------------------------
## 出力
##------------------------------------------------------------------------------

write.table(
  as.matrix(P_metropolis),
  file      = outfile1,
  sep       = "\t",
  row.names = FALSE,
  col.names = FALSE,
  quote     = FALSE
)

write.table(
  as.matrix(P_glauber),
  file      = outfile2,
  sep       = "\t",
  row.names = FALSE,
  col.names = FALSE,
  quote     = FALSE
)