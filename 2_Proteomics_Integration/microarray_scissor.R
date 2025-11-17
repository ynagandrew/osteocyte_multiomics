library(Scissor)
library(Seurat)
library(preprocessCore)
library(Matrix)
library(glue)
library(grid)
library(gridExtra)
library(qacBase)
library(rlang)

#Microarray Pre-processing
microarray <- read.delim('OSD-324/GLDS-324_microarray_2015-10-07_Pajevic_RMA-GENE-FULL_Group.txt', )
metadata <- read.delim2('OSD-324/s_GLDS-324.txt')

microarray$Gene.Symbol <- trimws(microarray$Gene.Symbol)
colnames(microarray)[2:13] <- metadata$Sample.Name

microarray <- microarray[!microarray$Gene.Symbol == "---",]
microarray <- microarray[!is.na(microarray$Gene.Symbol),]
microarray <- microarray[!duplicated(microarray$Gene.Symbol),]

rownames(microarray) <- microarray$Gene.Symbol

microarray <- microarray[2:13]
microarray <- as.matrix(microarray)

gravity <- integer(length(colnames(microarray)))
gravity[metadata$Factor.Value.Spaceflight. == "Space Flight"] <- 1
gravity <- as.numeric(gravity)

time <- as.factor(metadata$Factor.Value.Time)

time_2 <- integer(length(colnames(microarray)))
time_4 <- integer(length(colnames(microarray)))
time_6 <- integer(length(colnames(microarray)))

time_2[metadata$Factor.Value.Time == 2] <- 1
time_4[metadata$Factor.Value.Time == 4] <- 1
time_6[metadata$Factor.Value.Time == 6] <- 1

time_2 <- as.numeric(time_2)
time_4 <- as.numeric(time_4)
time_6 <- as.numeric(time_6)

tag_gravity <- c("Normal Gravity", "Microgravity")
# 
# Metacells Pre-processing
supercells <- readRDS("5_10_supercell.rds")
supercells <- supercells$seuratMC
supercells[["RNA"]] <- as(object = supercells[["RNA"]], Class = "Assay")

# 
# print("ligand")
# print(Sys.time())
# gravity_scissor <- Scissor(bulk_dataset = microarray,
#                  sc_dataset = supercells,
#                  phenotype = gravity,
#                  tag = tag_gravity,
#                  alpha = 0.3,
#                  family = "binomial",
#                  Save_file = "gravity.RData")
# saveRDS(gravity_scissor, file="gravity_scissor2.rds")

# gravity_scissor <- gravity_scissor2
# scissor <- gravity_scissor
# 
# Scissor_select <- rep(0, ncol(supercells))
# print(length(Scissor_select))
# names(Scissor_select) <- colnames(supercells)
# Scissor_select[scissor[["Scissor_pos"]]] <- 2
# Scissor_select[scissor[["Scissor_neg"]]] <- 1
# supercells <- AddMetaData(supercells, metadata = Scissor_select, col.name = "gravity2_scissor")


scissor_downstream <- function(scissor, poslabel, neglabel, filename, supercells) {
  
  
  pdf(file = glue("scissor2_plots/{filename}.pdf"), width=13, height=9)
  sc <- DimPlot(supercells, reduction = 'umap', group.by = 'gravity2_scissor', cols = c('grey','indianred1','royalblue'), pt.size = 1.2, order = c(2,1))
  grid.arrange(sc)
  dev.off()
  
  # load(glue("{filename}.RData"))
  #
  # numbers <- length(scissor$Scissor_pos) + length(scissor$Scissor_neg)
  # result2 <- reliability.test(X, Y, network, alpha = 0.2, family = "binomial",
  #                             cell_num = numbers, n = 10, nfold = 10)
  evaluate_summary <- evaluate.cell(glue("{filename}.RData"), scissor, FDR = 0.05, bootstrap_n = 100)
  
  return(c(supercells, evaluate_summary))
}

gravity_scissor <- readRDS("gravity_scissor2.rds")

scissor <- gravity_scissor

Scissor_select <- rep(0, ncol(supercells))
print(length(Scissor_select))
names(Scissor_select) <- colnames(supercells)
Scissor_select[scissor[["Scissor_pos"]]] <- 2
Scissor_select[scissor[["Scissor_neg"]]] <- 1
supercells <- AddMetaData(supercells, metadata = Scissor_select, col.name = "gravity2_scissor")

pdf(file = glue("scissor_plots/gravity.pdf"), width=9, height=9)
sc <- DimPlot(supercells, reduction = 'umap', group.by = "gravity2_scissor", cols = c('grey','indianred1','royalblue'), pt.size = 1.2, order = c(2,1))
grid.arrange(sc)
dev.off()

# gravity_scissor_results <- scissor_downstream(gravity_scissor, 1, 2, "gravity2_scissor", supercells)
