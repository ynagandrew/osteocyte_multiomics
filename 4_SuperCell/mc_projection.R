library(patchwork)
library(ggplot2)

mc_projection <- function(sc.obj,
                          mc.obj = NULL,
                          cell.membership = NULL,
                          metacell.label = NULL,
                          sc.label = NULL,
                          dims = c(1, 2),
                          sc.reduction = "umap",
                          mc.color = NULL,
                          sc.color = NULL,
                          alpha = 1,
                          pt_size = 1,
                          metric = "size",
                          continuous_metric = F) {
  
  if(is.null(mc.obj) & is.null(cell.membership)){
    stop("A membership vector should be provided either through the membership parameter or should be available in mc.obj as described in the documentation.")
  } 
  if(is.null(cell.membership)){
    membership <- mc.obj$SC$membership
  } else {
    membership <- cell.membership
  }  
  
  gamma <- mc.obj$SC$gamma
  knn <- mc.obj$SC$k.knn
  num_cells <- length(mc.obj$SC$graph.supercells)
  
  sc.obj$Metacell <- membership
  mc.obj <- mc.obj$seuratMC
  
  if(assertthat::is.string(sc.reduction)){
    # if sc_reduction does not exist compute pca and run UMAP:
    if(is.null(sc.obj@reductions[[sc.reduction]])){
      message("Low dimensionnal embessing not found in sc.obj")
      message("Computing PCA ...")
      sc.obj <- Seurat::NormalizeData(sc.obj, normalization.method = "LogNormalize")
      sc.obj <- Seurat::FindVariableFeatures(sc.obj, nfeatures = 2000)
      sc.obj <- Seurat::ScaleData(sc.obj)
      sc.obj <- Seurat::RunPCA(sc.obj, verbose = F)
      message("Running UMAP ...")
      sc.obj <- Seurat::RunUMAP(sc.obj, reduction = "pca", dims = c(1:30), n.neighbors = 15, verbose = F)
      scCoord <- Seurat::Embeddings(sc.obj@reductions[["umap"]])
    } else{
      scCoord <- Seurat::Embeddings(sc.obj@reductions[[sc.reduction]])
    } 
  } else if(!(is.data.frame(sc.reduction) | is.matrix(sc.reduction)) ){
    stop("sc.reduction should be a string indicating the name of the embedding to use in the reduction slot of sc.obj or a dataframe (or matrix) containing the components (columns) of single-cell embedding")
  } else{
    scCoord <- sc.reduction
  }  
  
  scCoordMetacell <-  cbind(scCoord, membership)
  
  centroids <- stats::aggregate(scCoord~membership, scCoord, mean) #should be taken from object slot
  rownames(centroids) <- centroids[,1]
  centroids <- centroids[colnames(mc.obj),]
  
  # if metacell.label and metric provided, mc.obj is mandatory if mc.obj not found membership mandatory and metric set to size which we compute from membership 
  # just add a message to warn the user 
  centroids[[metric]] <- mc.obj[[metric]][,1]
  
  if(is.null(metacell.label)) {
    metacell.label <- "MC"
    centroids[[metacell.label]] <- rep("red", length(mc.obj[[metric]]))
  } else {
    centroids[[metacell.label]] <- as.factor(mc.obj[[metacell.label]][,1])
    
    if(!is.null(sc.label)){
      if(metacell.label == sc.label){
        sc.obj@meta.data[, metacell.label] <- as.factor(sc.obj@meta.data[, metacell.label])
        centroids[[metacell.label]] <- factor(centroids[[metacell.label]], levels = levels(sc.obj@meta.data[, metacell.label]))
      } 
    } 
  }
  
  if(!is.null(sc.label)){
    
    scCoord <- data.frame(scCoord)
    scCoord[[sc.label]] <- sc.obj[[sc.label]][,1]
    
    p <- ggplot2::ggplot(scCoord,
                         ggplot2::aes_string(colnames(scCoord)[dims[1]],
                                             colnames(scCoord)[dims[2]],
                                             color = sc.label)) +
      ggplot2::geom_point(size=pt_size, alpha = alpha) +
      ggplot2::guides(color = ggplot2::guide_legend(override.aes = list(size = 2)))
    
    
  }else{
    p <- ggplot2::ggplot(data.frame(scCoord),
                         ggplot2::aes_string(colnames(scCoord)[dims[1]],
                                             colnames(scCoord)[dims[2]])) +
      ggplot2::geom_point(size=pt_size, color = "grey", alpha = 1)
  } 
  
  if(!continuous_metric){
    p <- p + ggplot2::geom_point(data=centroids, ggplot2::aes_string(colnames(centroids)[1 + dims[1]],
                                                                     colnames(centroids)[1 + dims[2]],
                                                                     fill = metacell.label, size = metric), colour="black", pch=21) 
  }else{
    p <- p + ggplot2::geom_point(data=centroids, ggplot2::aes_string(colnames(centroids)[1 + dims[1]],
                                                                     colnames(centroids)[1 + dims[2]],
                                                                     fill = metric), colour="black", pch=21, size=2) 
  } 
  
  if(!is.null(metacell.label) & !is.null(sc.label) & !is.null(mc.color)){
    if(metacell.label == sc.label & !is.null(mc.color)){
      sc.color = mc.color
    } 
  } 
  
  if(!is.null(metacell.label) & !is.null(sc.label) & !is.null(sc.color)){
    if(metacell.label == sc.label & !is.null(sc.color)){
      mc.color = sc.color
    } 
  } 
  
  if (!is.null(mc.color) & !continuous_metric) {
    p <- p + ggplot2::scale_fill_manual(values = mc.color) +  ggplot2::theme_classic() 
  }
  if (!is.null(sc.color)) {
    p <- p + ggplot2::scale_color_manual(values = sc.color) +  ggplot2::theme_classic()
  }
  
  p <- p + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                              panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
  

  # p <- p + ggplot2::guide_custom(title=glue("{gamma}"))
  
  return(p)
}