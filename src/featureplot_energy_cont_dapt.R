source("src/Functions.R")

args <- commandArgs(trailingOnly = TRUE)
sce_file    <- args[1]  # integrated SCE (with hexbin)
cont_file   <- args[2]  # cont seurat_landscaper.RData
dapt_file   <- args[3]  # DAPT seurat_landscaper.RData
outfile     <- args[4]

# Load integrated SCE (has hexbin structure)
load(sce_file)

# Load cont/DAPT Seurat objects into separate environments
cont_env <- new.env()
load(cont_file, envir = cont_env)
cont_seurat <- cont_env$seurat.integrated

dapt_env <- new.env()
load(dapt_file, envir = dapt_env)
dapt_seurat <- dapt_env$seurat.integrated

# Compute cont_DAPT_energy: per-cell energy from sub-analyses (negated)
cont_DAPT_energy <- rep(0, ncol(sce))
names(cont_DAPT_energy) <- colnames(sce)

cont_cells <- intersect(colnames(sce), colnames(cont_seurat))
dapt_cells <- intersect(colnames(sce), colnames(dapt_seurat))

cont_DAPT_energy[cont_cells] <- -cont_seurat$energy[cont_cells]
cont_DAPT_energy[dapt_cells] <- -dapt_seurat$energy[dapt_cells]

sce$cont_DAPT_energy <- cont_DAPT_energy

# Assign condition labels
existing_condition <- colData(sce)$condition
if (is.null(existing_condition) || length(existing_condition) == 0) {
    sce$condition <- ifelse(colnames(sce) %in% cont_cells, "cont",
                     ifelse(colnames(sce) %in% dapt_cells, "DAPT", NA))
} else {
    sce$condition <- gsub("-.*", "", existing_condition)
}

rm(cont_env, dapt_env, cont_seurat, dapt_seurat)

# ===== Hexbin energy diff =====
hex_ids <- sce@metadata$hexbin[[1]]
hex_positions <- as.data.frame(sce@metadata$hexbin[[2]])
hex_positions <- rownames_to_column(hex_positions, var = "hexID")

df <- data.frame(
  hexID = hex_ids,
  energy = sce$cont_DAPT_energy,
  condition = sce$condition
)
df$hexID <- as.character(df$hexID)
unique_hex <- sort(unique(hex_ids))
df$hexID <- match(hex_ids, unique_hex)

median_wide <- df %>%
  group_by(hexID, condition) %>%
  summarize(med_energy = median(energy, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = condition, values_from = med_energy) %>%
  filter(!is.na(cont) & !is.na(DAPT)) %>%
  mutate(diff = DAPT - cont)

hex_diff <- merge(median_wide, hex_positions, by = "hexID", all = FALSE)

# ===== Plots =====
p1 <- ggplot(hex_diff, aes(x = x, y = y, fill = cont)) +
  geom_hex(stat = "identity", na.rm = TRUE) +
  scale_fill_viridis_c() +
  coord_fixed() +
  theme_void() +
  theme(plot.background = element_rect(fill = "white", color = NA),
        panel.background = element_rect(fill = "white", color = NA),
        legend.position = "right",
        plot.title = element_text(hjust = 0.5, size = 14)) +
  labs(title = "Control Energy (Median)", fill = "cont")

p2 <- ggplot(hex_diff, aes(x = x, y = y, fill = DAPT)) +
  geom_hex(stat = "identity", na.rm = TRUE) +
  scale_fill_viridis_c() +
  coord_fixed() +
  theme_void() +
  theme(plot.background = element_rect(fill = "white", color = NA),
        panel.background = element_rect(fill = "white", color = NA),
        legend.position = "right",
        plot.title = element_text(hjust = 0.5, size = 14)) +
  labs(title = "DAPT Energy (Median)", fill = "DAPT")

p3 <- ggplot(hex_diff, aes(x = x, y = y, fill = diff)) +
  geom_hex(stat = "identity", na.rm = TRUE) +
  scale_fill_gradient2(low = "blue", mid = "gray80", high = "red", midpoint = 0) +
  coord_fixed() +
  theme_void() +
  theme(plot.background = element_rect(fill = "white", color = NA),
        panel.background = element_rect(fill = "white", color = NA),
        legend.position = "right",
        plot.title = element_text(hjust = 0.5, size = 14)) +
  labs(title = "Diff (DAPT - cont)", fill = "DAPT - cont")

png(file=outfile, width=1800, height=600)
grid.arrange(p1, p2, p3, ncol = 3)
dev.off()
