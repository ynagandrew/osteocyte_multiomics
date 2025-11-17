library(Seurat)
#library(glue)
#library(ggplot2)
#library(ggpubr)
#library(patchwork)
#library(VennDiagram)
library(dplyr)
library(SuperCell)
library(igraph)
library(glue)

for (knn in c(5, 15, 35)) {
  for (gamma in c(20, 50, 100)) {
    supercell <- readRDS(file=glue("test/{knn}_{gamma}_supercell.rds"))
    print(glue("knn: {knn} gamma: {gamma} num. metacells: {length(supercell$seuratMC@meta.data$size)}"))
    print(glue("size: avg {mean(supercell$seuratMC@meta.data$size)} sd {sd(supercell$seuratMC@meta.data$size)}"))
    print(glue("purity: avg {mean(supercell$seuratMC@meta.data$purity)} sd {sd(supercell$seuratMC@meta.data$purity)}"))
    print("")
  }
}