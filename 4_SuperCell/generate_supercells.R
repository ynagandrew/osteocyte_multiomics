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
print("hi, testing!")

source("/scratch/user/s4704475/SuperCell/mc_projection.R")

#allost <- readRDS("Osteocyte_v2025.rds")
# umap <- DimPlot(object = allost, reduction = "umap", group.by = 'ligand')
# pdf(file=glue("plots/5_umap.pdf"))
# grid.draw(umap)
# dev.off()
# 
# vln <- VlnPlot(allost, features = "Fyn")
# pdf(file=glue("plots/vln.pdf"))
# grid.draw(vln)
# dev.off()

# function from GfellerLab SIB Workshop w/ slight modifications
makeSeuratMC <- function(pbmc,
                         pc,
                         knn,
                         gamma,
                         genes = NULL,
                         metaFields = c("nCount_RNA", "nFeature_RNA", "time", "sample", "ligand", "stiffness", "percent.mt", "nCount_SCT", "nFeature_SCT"),
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
  
  # Separation, Compactness, INV adapted from Gfeller Lab Metacell Analysis Toolkit Source Code
  sc.reduction <- Seurat::Embeddings(allost, reduction = 'pca')
  cell.membership <- MC$membership
  
  # Separation
  print(Sys.time())
  print("Separation")
  print("-----")
  centroids <- stats::aggregate(x=sc.reduction, by= list(metacell = cell.membership), FUN = mean)
  dist_matrix <- as.matrix(dist(centroids[-1]))
  rownames(dist_matrix) <- centroids[,1]
  separation_distances <- apply(dist_matrix, 1, function(x) {x[order(x)[2]]})
  MC$separation <- separation_distances
  
  # Compactness
  print(Sys.time())
  print("Compactness")
  print("-----")
  centroids <- stats::aggregate(x=sc.reduction, by = list(metacell = cell.membership), FUN = var)
  compactness <- apply(centroids[, -1], 1, mean)
  names(compactness) <- centroids[,1]
  MC$compactness <- compactness
  
  # # Inner Normalized Variance
  # group.label = "membership"
  # assay = "RNA"
  # slot = "counts"
  # do.norm = T
  # NA_val = 1
  # scaling_factor = NULL
  # 
  # get_INV_val <- function(x, norm) {
  #   if(norm){
  #     median_val <- median(Matrix::rowSums(x))
  #     x_normalized <- (x / Matrix::rowSums(x)) * ifelse(is.null(scaling_factor), median_val, scaling_factor)
  #     result <- (1 / Matrix::colMeans(x_normalized)) * sparseMatrixStats::colVars(x_normalized)
  #   }else{
  #     result <- (1 / Matrix::colMeans(x)) * sparseMatrixStats::colVars(x)
  #   }
  #   return(result)
  # }
  # 
  # print(Sys.time())
  # print("INV")
  # print("-----")
  # data_matrix <- Matrix::t(Seurat::GetAssayData(pbmc, assay = assay, slot = slot))
  # result_list <- tapply(1:nrow(data_matrix), cell.membership, function(indices) {
  #   x <- data_matrix[indices, ]
  #   get_INV_val(x, norm = do.norm)
  # })
  # INV_val <- do.call(rbind, result_list)
  # INV_val[is.na(INV_val)] <- NA_val
  # 
  # # Compute the 95th percentile (quantile)
  # INV_val_qt <- apply(INV_val, 1, function(x) quantile(x, 0.95, na.rm = TRUE))
  # MC$INV <- INV_val_qt
  
  seuratMC <- supercell_2_Seurat(SC.GE = GE,MC,fields = c(metaFields,"purity","separation", "compactness"))
  seuratMC <- RunUMAP(seuratMC,dims = c(1:20))

  res <- seuratMC
  
  if (returnMC) {
    res <- list(seuratMC = seuratMC,SC = MC)
  }
  
  
  return(res)
}

knn_values <- c(5, 10, 15, 20, 35, 50, 100)
gamma_values <- c(10, 20, 50, 100, 250)



heatmappy <- function(input, direction, filename) {
  print(filename)
  test <- expand.grid(Knn=as.factor(knn_values), Gamma=as.factor(gamma_values))
  test$Z <- expand.grid(t(input))
  test$Z <- test$Z$Var1
  plot2 <- ggplot(test, aes(x=Gamma, y=Knn, fill = Z)) +
    geom_tile() + geom_text(aes(label=round(Z, 3)), color="white") +
    scale_fill_viridis_c(, direction = direction) + coord_fixed() + 
    guides(fill = guide_colourbar(barwidth = 0.5,
                                  barheight = 10, title=filename))
  
  ggsave(glue("{filename}_heatmap.png"), width = 5, height = 6)
}

lineplotty <- function(input, name) {
  test <- expand.grid(Knn=as.factor(knn_values), Gamma=as.factor(gamma_values))
  test$Z <- expand.grid(t(input))
  test$Z <- test$Z$Var1
  plot2 <- ggplot(test, aes(x=Gamma, y=Z, group=Knn)) +
    geom_line(aes(color=Knn)) + 
    geom_point(aes(color=Knn)) + 
    labs(title=name, y=name)
  plot2
}

# purity_grid <- matrix(1:35, nrow = 5, dimnames = list(gamma_values, knn_values))
# separation_grid <- matrix(1:35, nrow = 5, dimnames = list(gamma_values, knn_values))
# compactness_grid <- matrix(1:35, nrow = 5, dimnames = list(gamma_values, knn_values))
# 
# purity_grid <- readRDS("purity_grid.rds")
# separation_grid <- readRDS("separation_grid.rds")
# compactness_grid <- readRDS("compactness_grid.rds")

density_grid <- matrix(1:35, nrow = 5, dimnames = list(gamma_values, knn_values))


for (knn in knn_values) {
  for (gamma in gamma_values) {
    
    print(Sys.time())
    print(knn)
    print(gamma)
    supercell <- readRDS(file=glue("supercells/{knn}_{gamma}_supercell.rds"))
    density_grid[as.character(gamma), as.character(knn)] <- mean(supercell$seuratMC@meta.data[["nFeature_RNA"]])/21071

    # # Generate Metacells w/ different KNN values.
    # supercell <- makeSeuratMC(allost, 20, knn, gamma, returnMC=T)
    # supercell$seuratMC <- RunUMAP(supercell$seuratMC,dims = c(1:20))
    # 
    # saveRDS(supercell, file = glue("supercells/{knn}_{gamma}_supercell.rds"))
    # 
    # plot0 <- mc_projection(allost, supercell, metacell.label = "sample", sc.label="sample", pt_size=0.01)
    # png(file=glue("plot_metacell/{knn}_{gamma}_supercell_UMAP.png"), width=1080, height=1080)
    # grid.draw(plot0)
    # dev.off()
    
    ## Iterate through Metacells w/ different KNN values & plot QC stats
    # purity_df <- supercell$seuratMC@meta.data[c("sample", "purity")]
    # purity_per_sample <- aggregate(purity_df$purity, list(purity_df$sample), mean)
    # purity_score <- mean(purity_per_sample$x)
    # purity_grid[as.character(gamma), as.character(knn)] <- purity_score
    # 
    # separation_df <- supercell$seuratMC@meta.data[c("sample", "separation")]
    # separation_per_sample <- aggregate(separation_df$separation, list(separation_df$sample), mean)
    # separation_score <- mean(separation_per_sample$x)
    # separation_grid[as.character(gamma), as.character(knn)] <- separation_score
    # 
    # compactness_df <- supercell$seuratMC@meta.data[c("sample", "compactness")]
    # 
    # compactness_df <- compactness_df[complete.cases(compactness_df),]
    # compactness_per_sample <- aggregate(compactness_df$compactness, list(compactness_df$sample), mean)
    # compactness_score <- mean(compactness_per_sample$x)
    # compactness_grid[as.character(gamma), as.character(knn)] <- compactness_score
    # 
    # saveRDS(purity_grid, "purity_grid.rds")
    # saveRDS(compactness_grid, "compactness_grid.rds")
    # saveRDS(separation_grid, "separation_grid.rds")
    }
}

# heatmappy(separation_grid, "separation.pdf")
# heatmappy(purity_grid, "purity.pdf")
# heatmappy(compactness_grid, "compactness")
heatmappy(density_grid, 1,"density")
