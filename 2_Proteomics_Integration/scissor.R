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

# time_1 <- integer(length(colnames(proteomics)))
# time_3 <- integer(length(colnames(proteomics)))
# time_7 <- integer(length(colnames(proteomics)))
# 
# ligand <- integer(length(colnames(proteomics)))
# gravity <- integer(length(colnames(proteomics)))
# 
# 
# time_1[grep("Day.1", colnames(proteomics))] <- 1
# time_3[grep("Day.3", colnames(proteomics))] <- 1
# time_7[grep("Day.7", colnames(proteomics))] <- 1
# 
# ligand[grep("with", colnames(proteomics))] <- 1
# 
# gravity[grep("Normal", colnames(proteomics))] <- 1
# 
# 
# time_1 <- as.numeric(as.character(time_1))
# time_3 <- as.numeric(as.character(time_3))
# time_7 <- as.numeric(as.character(time_7))
# 
# ligand <- as.numeric(as.character(ligand))
# gravity <- as.numeric(as.character(gravity))
# 
# 
# tag_1 <- c("Not Day 1", "Day 1")
# tag_3 <- c("Not Day 3", "Day 3")
# tag_7 <- c("Not Day 7", "Day 7")
# 
# tag_ligand = c("No Ligand", "With Ligand")
# tag_gravity = c("Microgravity", "Normal Gravity")
# supercells[["RNA"]] <- as(object = supercells[["RNA"]], Class = "Assay")
# proteomics <- as.matrix(proteomics)
# 
# print("ligand")
# print(Sys.time())
# ligand_scissor <- Scissor(bulk_dataset = proteomics,
#                  sc_dataset = supercells,
#                  phenotype = ligand,
#                  tag = tag_ligand,
#                  alpha = 0.5,
#                  family = "binomial",
#                  Save_file = "ligand.RData")
# saveRDS(ligand_scissor, file="ligand_scissor.rds")
# 
# 
# print("gravity")
# print(Sys.time())
# gravity_scissor <- Scissor(bulk_dataset = proteomics,
#                   sc_dataset = supercells,
#                   phenotype = gravity,
#                   tag = tag_gravity,
#                   alpha = 0.5,
#                   family = "binomial")
# saveRDS(gravity_scissor, file="gravity_scissor.rds")
# 
# print("time_1")
# print(Sys.time())
# time1_scissor <- Scissor(bulk_dataset = proteomics,
#                            sc_dataset = supercells,
#                            phenotype = time_1,
#                            tag = tag_1,
#                            alpha = 0.5,
#                            family = "binomial")
# # saveRDS(time1_scissor, file="time_1_scissor.rds")
# 
# print("time_3")
# print(Sys.time())
# time3_scissor <- Scissor(bulk_dataset = proteomics,
#                          sc_dataset = supercells,
#                          phenotype = time_3,
#                          tag = tag_3,
#                          alpha = 0.5,
#                          family = "binomial")
# saveRDS(time3_scissor, file="time_3_scissor.rds")
# 
# print("time_7")
# print(Sys.time())
# time7_scissor <- Scissor(bulk_dataset = proteomics,
#                          sc_dataset = supercells,
#                          phenotype = time_7,
#                          tag = tag_7,
#                          alpha = 0.5,
#                          family = "binomial")
# saveRDS(time7_scissor, file="time_7_scissor.rds")

# all_scissors = c("ligand_scissor" = ligand_scissor, "gravity_scissor" = gravity_scissor,
#                  "time_1_scissor" = time_1_scissor, "time3_scissor" = time_3_scissor,
#                  "time_7_scissor" = time_7_scissor)
# 
# for (scissor in all_scissors) {
#   print(deparse(substitute(scissor)))
#   Scissor_select <- rep(0, ncol(supercells))
#   print(length(Scissor_select))
#   names(Scissor_select) <- colnames(supercells)
#   Scissor_select[scissor[["Scissor_pos"]]] <- 2
#   Scissor_select[scissor[["Scissor_neg"]]] <- 1
#   supercells <- AddMetaData(supercells, metadata = Scissor_select, col.name = deparse(substitute(scissor)))
# }


scissor = time_7_scissor

Scissor_select <- rep(0, ncol(supercells))
print(length(Scissor_select))
names(Scissor_select) <- colnames(supercells)
Scissor_select[scissor[["Scissor_pos"]]] <- 2
Scissor_select[scissor[["Scissor_neg"]]] <- 1
supercells <- AddMetaData(supercells, metadata = Scissor_select, col.name = "time_7_scissor")

pdf(file = glue("scissor_plots/time7.pdf"), width=9, height=9)
sc <- DimPlot(supercells, reduction = 'umap', group.by = "time_7_scissor", cols = c('grey','indianred1','royalblue'), pt.size = 1.2, order = c(2,1))
grid.arrange(sc)
dev.off()


# all_names = c("ligand_scissor", "gravity_scissor",
#                  "time_1_scissor", "time_3_scissor",
#                  "time_7_scissor")
# 
# for (names in all_names) {
#   pdf(file = glue("scissor_plots/{names}.pdf"), width=13, height=9)
#   sc <- DimPlot(supercells, reduction = 'umap', group.by = names, cols = c('grey','indianred1','royalblue'), pt.size = 1.2, order = c(2,1))
#   grid.arrange(sc)
#   dev.off()
# }

# scissor_downstream <- function(scissor, poslabel, neglabel, filename, supercells) {
# 
# 
#   pdf(file = glue("scissor_plots/{filename}.pdf"), width=13, height=9)
#   sc <- DimPlot(supercells, reduction = 'umap', group.by = 'scissor', cols = c('grey','indianred1','royalblue'), pt.size = 1.2, order = c(2,1))
#   grid.arrange(sc)
#   dev.off()
# 
#   # load(glue("{filename}.RData"))
#   #
#   # numbers <- length(scissor$Scissor_pos) + length(scissor$Scissor_neg)
#   # result2 <- reliability.test(X, Y, network, alpha = 0.2, family = "binomial",
#   #                             cell_num = numbers, n = 10, nfold = 10)
#   evaluate_summary <- evaluate.cell(glue("{filename}.RData"), scissor, FDR = 0.05, bootstrap_n = 100)
# 
#   return(c(supercells, evaluate_summary))
# }

# ligand_scissor <- readRDS("ligand_scissor.rds")
# 
# ligand_scissor_results <- scissor_downstream(ligand_scissor, 1, 2, "ligand", supercells)
# saveRDS(ligand_scissor_results, "ligand_scissor_results.rds")
