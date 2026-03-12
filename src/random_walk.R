source("src/Functions.R")

# Parameter
infile1  <- commandArgs(trailingOnly = TRUE)[1]
infile2  <- commandArgs(trailingOnly = TRUE)[2]
outfile1 <- commandArgs(trailingOnly = TRUE)[3]
outfile2 <- commandArgs(trailingOnly = TRUE)[4]
# infile1 <- 'plot/hpbase/integrated/P_metropolis.tsv'
# infile2 <- 'plot/hpbase/integrated/P_glauber.tsv'

# Settings
n_steps <- 10   # 1,2,5,10 など自由に
# tol / max_iter はもう不要（途中で止めない）

# Loading (row-stochastic P[i,j] = i -> j)
P_m <- as.matrix(read.table(infile1, header = FALSE))
P_g <- as.matrix(read.table(infile2, header = FALSE))

# -------------------------
# Random walk (metropolis)
# -------------------------

P_m_list <- vector("list", n_steps + 1)
P_m_list[[1]] <- diag(nrow(P_m))   # P^0 = I

P_m_cur <- P_m_list[[1]]
for (iter in seq_len(n_steps)) {
  P_m_cur <- P_m_cur %*% P_m
  P_m_list[[iter + 1]] <- P_m_cur
}
names(P_m_list) <- paste0("step_", 0:n_steps)

print(paste("Metropolis steps:", n_steps))

# -------------------------
# Random walk (glauber)
# -------------------------

P_g_list <- vector("list", n_steps + 1)
P_g_list[[1]] <- diag(nrow(P_g))   # P^0 = I

P_g_cur <- P_g_list[[1]]
for (iter in seq_len(n_steps)) {
  P_g_cur <- P_g_cur %*% P_g
  P_g_list[[iter + 1]] <- P_g_cur
}
names(P_g_list) <- paste0("step_", 0:n_steps)

print(paste("Glauber steps:", n_steps))

# -------------------------
# Save (list output)
# -------------------------

save(P_m_list, file = outfile1)
save(P_g_list, file = outfile2)
