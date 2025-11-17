library(Seurat)
library(glue)
library(ggplot2)
#library(ggpubr)
#library(patchwork)
#library(VennDiagram)
library(dplyr)
library(SuperCell)
library(stats)
library(grid)
library(colorspace)
library(viridis)

source("/scratch/user/s4704475/SuperCell/mc_projection.R")
source("/scratch/user/s4704475/SuperCell/generate_supercells.R")

pbmc <- readRDS(file = "Osteocyte_v2025.rds")
knn_values <- c(5, 10, 15, 20, 35, 50, 100)
gamma_values <- c(10, 20, 50, 100, 250)

purity_grid <- matrix(1:35, nrow = 5, dimnames = list(gamma_values, knn_values))
separation_grid <- matrix(1:35, nrow = 5, dimnames = list(gamma_values, knn_values))
compactness_grid <- matrix(1:35, nrow = 5, dimnames = list(gamma_values, knn_values))

time_purity_grid <- matrix(1:35, nrow = 5, dimnames = list(gamma_values, knn_values))
ligand_purity_grid <- matrix(1:35, nrow = 5, dimnames = list(gamma_values, knn_values))
stiffness_purity_grid <- matrix(1:35, nrow = 5, dimnames = list(gamma_values, knn_values))
name_purity_grid <- matrix(1:35, nrow = 5, dimnames = list(gamma_values, knn_values))

nCount_RNA_grid <- matrix(1:35, nrow = 5, dimnames = list(gamma_values, knn_values))

for (knn in knn_values) {
  for (gamma in gamma_values) {
    print(Sys.time())
    print(knn)
    print(gamma)
    
    supercell <- readRDS(file=glue("supercells/{knn}_{gamma}_supercell.rds"))
    
    supercell$seuratMC$time_purity <- supercell_purity(clusters = pbmc$time,
                                  supercell_membership = supercell$SC$membership)

    supercell$seuratMC$ligand_purity <- supercell_purity(clusters = pbmc$ligand,
                                  supercell_membership = supercell$SC$membership)

    supercell$seuratMC$stiffness_purity <- supercell_purity(clusters = pbmc$stiffness,
                                  supercell_membership = supercell$SC$membership)
    
    purity_df <- supercell$seuratMC@meta.data[c("sample", "purity")]
    purity_per_sample <- aggregate(purity_df$purity, list(purity_df$sample), mean)
    purity_score <- mean(purity_per_sample$x)
    purity_grid[as.character(gamma), as.character(knn)] <- purity_score
    
    separation_df <- supercell$seuratMC@meta.data[c("sample", "separation")]
    separation_per_sample <- aggregate(separation_df$separation, list(separation_df$sample), mean)
    separation_score <- mean(separation_per_sample$x)
    separation_grid[as.character(gamma), as.character(knn)] <- separation_score
    
    compactness_df <- supercell$seuratMC@meta.data[c("sample", "compactness")]
    
    compactness_df <- compactness_df[complete.cases(compactness_df),]
    compactness_per_sample <- aggregate(compactness_df$compactness, list(compactness_df$sample), mean)
    compactness_score <- mean(compactness_per_sample$x)
    compactness_grid[as.character(gamma), as.character(knn)] <- compactness_score

    time_purity_df <- supercell$seuratMC@meta.data[c("sample", "time_purity")]
    time_purity_per_sample <- aggregate(time_purity_df$time_purity, list(time_purity_df$sample), mean)
    time_purity_score <- mean(time_purity_per_sample$x)
    time_purity_grid[as.character(gamma), as.character(knn)] <- time_purity_score

    ligand_purity_df <- supercell$seuratMC@meta.data[c("sample", "ligand_purity")]
    ligand_purity_per_sample <- aggregate(ligand_purity_df$ligand_purity, list(ligand_purity_df$sample), mean)
    ligand_purity_score <- mean(ligand_purity_per_sample$x)
    ligand_purity_grid[as.character(gamma), as.character(knn)] <- ligand_purity_score

    stiffness_purity_df <- supercell$seuratMC@meta.data[c("sample", "stiffness_purity")]
    stiffness_purity_per_sample <- aggregate(stiffness_purity_df$stiffness_purity, list(stiffness_purity_df$sample), mean)
    stiffness_purity_score <- mean(stiffness_purity_per_sample$x)
    stiffness_purity_grid[as.character(gamma), as.character(knn)] <- stiffness_purity_score
    
    nCount_RNA_df <- as.numeric(supercell$SC$nCount_RNA)
    nCount_RNA_grid[as.character(gamma), as.character(knn)] <- mean(nCount_RNA_df)
    
    saveRDS(supercell, file = glue("supercells/{knn}_{gamma}_supercell.rds"))
    
  }
}

heatmappy(separation_grid, 1, "separation")
heatmappy(purity_grid, 1, "purity")
heatmappy(compactness_grid, -1, "compactness")

heatmappy(time_purity_grid, 1, "time_purity")
heatmappy(ligand_purity_grid, 1, "ligand_purity")
heatmappy(stiffness_purity_grid, 1, "stiffness_purity")

