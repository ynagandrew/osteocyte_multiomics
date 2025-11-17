library(Seurat)
library(ggplot2)
library(dplyr)

Osteocyte_v2025 <- readRDS("~/Desktop/Scratch/s4704475/SuperCell/Osteocyte_v2025.rds")
`5_20_supercell` <- readRDS("~/Desktop/Scratch/s4704475/SuperCell/supercells/5_20_supercell.rds")

# Extract the data matrices (log-normalized)
mat1 <- GetAssayData(Osteocyte_v2025, assay = "SCT", slot = "data")
mat2 <- GetAssayData(`5_20_supercell`[[1]], assay = "RNA", slot = "data")

# Flatten into long numeric vectors
expr1 <- as.vector(mat1)
expr2 <- as.vector(mat2)

# Combine into one data frame
df <- data.frame(
  expression = c(expr1, expr2),
  dataset = c(
    rep("Osteocyte_v2025", length(expr1)),
    rep("5_20_supercell", length(expr2))
  )
)

# Optionally filter out zeros if you only want nonzero expression
# df <- df %>% filter(expression > 0)

set.seed(123)
df_sample <- df %>% sample_n(min(100000, nrow(df)))

ggplot(df_sample, aes(x = dataset, y = expression, fill = dataset)) +
  geom_violin(trim = FALSE, scale = "width", width = 0.8) +
  geom_boxplot(width = 0.15, outlier.shape = NA, alpha = 0.4) +
  theme_classic(base_size = 14) +
  labs(y = "Normalized expression (log1p scale)",
       x = NULL,
       title = "Overall expression distribution (all genes × all cells)") +
  theme(legend.position = "none") +
  scale_y_continuous(trans = "log1p")  # optional, if not already log-normalized
