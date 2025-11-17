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

print("distinct_similar_analysis!")

source("/scratch/user/s4704475/SuperCell/mc_projection.R")
source("/scratch/user/s4704475/SuperCell/generate_supercells.R")
# knn_values <- c(5, 10, 15, 20, 35, 50, 100)
# gamma_values <- c(10, 20, 50, 100, 250)
# 
# purity_grid_5_S_LIN <- matrix(1:35, nrow = 5, dimnames = list(gamma_values, knn_values))
# separation_grid_5_S_LIN <- matrix(1:35, nrow = 5, dimnames = list(gamma_values, knn_values))
# compactness_grid_5_S_LIN <- matrix(1:35, nrow = 5, dimnames = list(gamma_values, knn_values))
# 
# purity_grid_72_S_LIN <- matrix(1:35, nrow = 5, dimnames = list(gamma_values, knn_values))
# separation_grid_72_S_LIN <- matrix(1:35, nrow = 5, dimnames = list(gamma_values, knn_values))
# compactness_grid_72_S_LIN <- matrix(1:35, nrow = 5, dimnames = list(gamma_values, knn_values))
# 
# purity_grid_5_L_CYC <- matrix(1:35, nrow = 5, dimnames = list(gamma_values, knn_values))
# separation_grid_5_L_CYC <- matrix(1:35, nrow = 5, dimnames = list(gamma_values, knn_values))
# compactness_grid_5_L_CYC <- matrix(1:35, nrow = 5, dimnames = list(gamma_values, knn_values))
# 
# purity_grid_72_L_CYC <- matrix(1:35, nrow = 5, dimnames = list(gamma_values, knn_values))
# separation_grid_72_L_CYC <- matrix(1:35, nrow = 5, dimnames = list(gamma_values, knn_values))
# compactness_grid_72_L_CYC <- matrix(1:35, nrow = 5, dimnames = list(gamma_values, knn_values))
# 
# for (knn in knn_values) {
#   for (gamma in gamma_values) {
#     
#     print(Sys.time())
#     print(knn)
#     print(gamma)
#     
#     supercell <- readRDS(file=glue("test/{knn}_{gamma}_supercell.rds"))
#     
#     purity_df <- supercell$seuratMC@meta.data[c("sample", "purity")]
#     purity_per_sample <- aggregate(purity_df$purity, list(purity_df$sample), mean)
# 
#     purity_5_S_LIN <- purity_per_sample$x[purity_per_sample$Group.1 == "SLIN_5hr"]
#     purity_grid_5_S_LIN[as.character(gamma), as.character(knn)] <- purity_5_S_LIN
# 
#     purity_72_S_LIN <- purity_per_sample$x[purity_per_sample$Group.1 == "SLIN_72hr"]
#     purity_grid_72_S_LIN[as.character(gamma), as.character(knn)] <- purity_72_S_LIN
# 
#     purity_5_L_CYC <- purity_per_sample$x[purity_per_sample$Group.1 == "LCYC_5hr"]
#     purity_grid_5_L_CYC[as.character(gamma), as.character(knn)] <- purity_5_L_CYC
# 
#     purity_72_L_CYC <- purity_per_sample$x[purity_per_sample$Group.1 == "LCYC_72hr"]
#     purity_grid_72_L_CYC[as.character(gamma), as.character(knn)] <- purity_72_L_CYC
#     
#     compactness_df <- supercell$seuratMC@meta.data[c("sample", "compactness")]
#     compactness_per_sample <- aggregate(compactness_df$compactness, list(compactness_df$sample), mean)
#     
#     compactness_5_S_LIN <- compactness_per_sample$x[compactness_per_sample$Group.1 == "SLIN_5hr"]
#     compactness_grid_5_S_LIN[as.character(gamma), as.character(knn)] <- compactness_5_S_LIN
#     
#     compactness_72_S_LIN <- compactness_per_sample$x[compactness_per_sample$Group.1 == "SLIN_72hr"]
#     compactness_grid_72_S_LIN[as.character(gamma), as.character(knn)] <- compactness_72_S_LIN
#     
#     compactness_5_L_CYC <- compactness_per_sample$x[compactness_per_sample$Group.1 == "LCYC_5hr"]
#     compactness_grid_5_L_CYC[as.character(gamma), as.character(knn)] <- compactness_5_L_CYC
#     
#     compactness_72_L_CYC <- compactness_per_sample$x[compactness_per_sample$Group.1 == "LCYC_72hr"]
#     compactness_grid_72_L_CYC[as.character(gamma), as.character(knn)] <- compactness_72_L_CYC
#     
#     separation_df <- supercell$seuratMC@meta.data[c("sample", "separation")]
#     separation_per_sample <- aggregate(separation_df$separation, list(separation_df$sample), mean)
#     
#     separation_5_S_LIN <- separation_per_sample$x[separation_per_sample$Group.1 == "SLIN_5hr"]
#     separation_grid_5_S_LIN[as.character(gamma), as.character(knn)] <- separation_5_S_LIN
#     
#     separation_72_S_LIN <- separation_per_sample$x[separation_per_sample$Group.1 == "SLIN_72hr"]
#     separation_grid_72_S_LIN[as.character(gamma), as.character(knn)] <- separation_72_S_LIN
#     
#     separation_5_L_CYC <- separation_per_sample$x[separation_per_sample$Group.1 == "LCYC_5hr"]
#     separation_grid_5_L_CYC[as.character(gamma), as.character(knn)] <- separation_5_L_CYC
#     
#     separation_72_L_CYC <- separation_per_sample$x[separation_per_sample$Group.1 == "LCYC_72hr"]
#     separation_grid_72_L_CYC[as.character(gamma), as.character(knn)] <- separation_72_L_CYC
#   }
# }
# 
# lineplotty(compactness_grid_5_S_LIN, "compactness_5_S_LIN")
# lineplotty(compactness_grid_5_L_CYC, "compactness_5_L_CYC")
# lineplotty(compactness_grid_72_S_LIN, "compactness_72_S_LIN")
# lineplotty(compactness_grid_72_L_CYC, "compactness_72_L_CYC")
# 
# saveRDS(purity_grid_5_L_CYC, "distinct_sample_metacell_test/purity_5_L_CYC.rds")
# saveRDS(purity_grid_5_S_LIN, "distinct_sample_metacell_test/purity_5_S_LIN.rds")
# saveRDS(purity_grid_72_L_CYC, "distinct_sample_metacell_test/purity_72_L_CYC.rds")
# saveRDS(purity_grid_72_S_LIN, "distinct_sample_metacell_test/purity_72_S_LIN.rds")
# #
# saveRDS(separation_grid_5_L_CYC, "distinct_sample_metacell_test/separation_5_L_CYC.rds")
# saveRDS(separation_grid_5_S_LIN, "distinct_sample_metacell_test/separation_5_S_LIN.rds")
# saveRDS(separation_grid_72_L_CYC, "distinct_sample_metacell_test/separation_72_L_CYC.rds")
# saveRDS(separation_grid_72_S_LIN, "distinct_sample_metacell_test/separation_72_S_LIN.rds")
# #
# saveRDS(compactness_grid_5_L_CYC, "distinct_sample_metacell_test/compactness_5_L_CYC.rds")
# saveRDS(compactness_grid_5_S_LIN, "distinct_sample_metacell_test/compactness_5_S_LIN.rds")
# saveRDS(compactness_grid_72_L_CYC, "distinct_sample_metacell_test/compactness_72_L_CYC.rds")
# saveRDS(compactness_grid_72_S_LIN, "distinct_sample_metacell_test/compactness_72_S_LIN.rds")

purity_5_L_CYC <- readRDS(file = "distinct_sample_metacell_test/purity_5_L_CYC.rds")
purity_72_L_CYC <- readRDS(file = "distinct_sample_metacell_test/purity_72_L_CYC.rds")
purity_5_S_LIN <- readRDS(file = "distinct_sample_metacell_test/purity_5_S_LIN.rds")

different_purity <- purity_5_L_CYC - purity_72_L_CYC
similar_purity <- purity_5_L_CYC - purity_5_S_LIN

lineplotty(different_purity, name = "5_L_CYC - 72_L_CYC")
lineplotty(similar_purity, name = "5_L_CYC - 5_S_LIN")

compactness_5_L_CYC <- readRDS(file = "distinct_sample_metacell_test/compactness_5_L_CYC.rds")
compactness_72_L_CYC <- readRDS(file = "distinct_sample_metacell_test/compactness_72_L_CYC.rds")
compactness_5_S_LIN <- readRDS(file = "distinct_sample_metacell_test/compactness_5_S_LIN.rds")

different_compactness <- compactness_5_L_CYC - compactness_72_L_CYC
similar_compactness <- compactness_5_L_CYC - compactness_5_S_LIN

lineplotty(different_compactness, name = "compactness 5_L_CYC - 72_L_CYC")
lineplotty(similar_compactness, name = "compactness 5_L_CYC - 5_S_LIN")

separation_5_L_CYC <- readRDS(file = "distinct_sample_metacell_test/separation_5_L_CYC.rds")
separation_72_L_CYC <- readRDS(file = "distinct_sample_metacell_test/separation_72_L_CYC.rds")
separation_5_S_LIN <- readRDS(file = "distinct_sample_metacell_test/separation_5_S_LIN.rds")

different_separation <- separation_5_L_CYC - separation_72_L_CYC
similar_separation <- separation_5_L_CYC - separation_5_S_LIN

lineplotty(different_separation, name = "separation 5_L_CYC - 72_L_CYC")
lineplotty(similar_separation, name = "separation 5_L_CYC - 5_S_LIN")