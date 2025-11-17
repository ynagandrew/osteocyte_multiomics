library(Seurat)
library(ggplot2)
library(grid)
library(gridExtra)
library(MOFA2)
library(data.table)
library(tidyr)
library(dplyr)
library(reticulate)

proteomics <- readRDS("proteomics_reduct.rds")
supercells <- readRDS("../SuperCell/supercells/5_10_supercell.rds")
supercells <- supercells$seuratMC

pseudo_supercells <- AggregateExpression(supercells, assays = "RNA", return.seurat = T, group.by = "sample")