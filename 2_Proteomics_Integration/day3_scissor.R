library(Scissor)
library(Seurat)
library(preprocessCore)
library(Matrix)
library(glue)
library(grid)
library(gridExtra)
library(qacBase)
library(rlang)

proteomics <- readRDS("../day1_3_7.rds")
supercells <- readRDS("5_10_supercell.rds")
supercells <- supercells$seuratMC

proteomics <- proteomics[grep("Day.3", colnames(proteomics))]

ligand <- integer(length(colnames(proteomics)))
gravity <- integer(length(colnames(proteomics)))

ligand[grep("with", colnames(proteomics))] <- 1

gravity[grep("Normal", colnames(proteomics))] <- 1

ligand <- as.numeric(as.character(ligand))
gravity <- as.numeric(as.character(gravity))

tag_ligand = c("No Ligand", "With Ligand")
tag_gravity = c("Microgravity", "Normal Gravity")
supercells[["RNA"]] <- as(object = supercells[["RNA"]], Class = "Assay")
proteomics <- as.matrix(proteomics)

print("gravity")
print(Sys.time())
gravity_scissor <- Scissor(bulk_dataset = proteomics,
                  sc_dataset = supercells,
                  phenotype = gravity,
                  tag = tag_gravity,
                  alpha = 0.01,
                  family = "binomial")

alphas = c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0)

for (alphaa in alphas) {
  print("ligand")
  print(alphaa)
  print(Sys.time())
  ligand_scissor <- Scissor(bulk_dataset = proteomics,
                            sc_dataset = supercells,
                            phenotype = ligand,
                            tag = tag_ligand,
                            alpha = alphaa,
                            family = "binomial")
  
  Scissor_select <- rep(0, ncol(supercells))
  print(length(Scissor_select))
  names(Scissor_select) <- colnames(supercells)
  Scissor_select[ligand_scissor[["Scissor_pos"]]] <- 2
  Scissor_select[ligand_scissor[["Scissor_neg"]]] <- 1
  supercells <- AddMetaData(supercells, metadata = Scissor_select, col.name = "ligand_scissor")
  
  pdf(file = glue("day3_scissor_plots/ligand_scissor_{alphaa}.pdf"), width=13, height=9)
  sc <- DimPlot(supercells, reduction = 'umap', group.by = "ligand_scissor", cols = c('grey','indianred1','royalblue'), pt.size = 1.2, order = c(2,1))
  grid.arrange(sc)
  dev.off()
}
