source("src/Functions.R")

# Parameter
args <- commandArgs(trailingOnly = TRUE)
input_rdata <- args[1]
outfile1    <- args[2]
outfile2    <- args[3]

##------------------------------------------------------------------------------
## ロード
##------------------------------------------------------------------------------

load(input_rdata)

E         <- seurat.integrated@misc$landscaper$E
Allstates <- seurat.integrated@misc$landscaper$Allstates
cell_time <- seurat.integrated$cell_time
state_idx <- seurat.integrated$state_idx

P    <- nrow(Allstates)

##------------------------------------------------------------------------------
## 状態ごとの平均時間 time_state（長さ P）を作成
##------------------------------------------------------------------------------

time_state <- rep(NA_real_, P)

tmp <- tapply(cell_time, state_idx, mean)  # names(tmp) は state_idx の値（1..P）

idx_states_obs <- as.integer(names(tmp))
time_state[idx_states_obs] <- as.numeric(tmp)

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
  time_state = time_state,
  tol        = tol_h
)

P_glauber <- transition_matrix_from_E(
  E          = E,
  S          = Allstates,
  beta       = 1,
  kernel     = "glauber",
  time_state = time_state,
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