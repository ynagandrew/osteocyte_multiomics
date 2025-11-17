library(Seurat)
library(SuperCell)
allost <- readRDS('Ost_all.rds')

ost_list <- SplitObject(allost, split.by = "time")

makeSeuratMC <- function(pbmc,
                         pc,
                         knn,
                         gamma,
                         genes = NULL,
                         metaFields = c("nCount_RNA", "nFeature_RNA", "time", "sample", "ligand", "stiffness", "percent.mt", "nCount_SCT", "nFeature_SCT", "name"),
                         returnMC = F) {
  
  if(is.null(genes)) {
    genes <- VariableFeatures(pbmc)
  }
  
  MC <- SCimplify(GetAssayData(pbmc,layer = "data"),  # normalized gene expression matrix 
                  n.pc = pc,
                  k.knn = knn, # number of nearest neighbors to build kNN network
                  gamma = gamma, # graining level
                  genes.use = genes )# will be the ones used for integration if input is seurat integrated data
  
  
  MC$purity <- supercell_purity(clusters = pbmc$sample,
                                supercell_membership = MC$membership)
  
  
  for (m in metaFields) {
    MC[[m]]<- supercell_assign(clusters = pbmc@meta.data[,m], # single-cell assigment to cell lines (clusters)
                               supercell_membership = MC$membership, # single-cell assignment to super-cells
                               method = "absolute")
  }
  
  GE <- supercell_GE(as.matrix(GetAssayData(pbmc,layer = "data")),groups = MC$membership)
  
  seuratMC <- supercell_2_Seurat(SC.GE = GE,MC,fields = c(metaFields,"purity"))
  
  seuratMC <- RunUMAP(seuratMC,dims = c(1:20))
  
  res <- seuratMC
  
  if (returnMC) {
    res <- list(seuratMC = seuratMC,SC = MC)
  }
  
  
  return(res)
}

hr_5 <- makeSeuratMC(ost_list[[1]], 20, 5, 20, returnMC = TRUE)
hr_72 <- makeSeuratMC(ost_list[[2]], 20, 5, 20, returnMC = TRUE)