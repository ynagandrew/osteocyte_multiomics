library(Seurat)
library(CIDER)

# seurat_integrate = readRDS("seurat.combined.rds")
# seurat_integrate <- seurat_integrate[[1]]
# seu.integrated <- hdbscan.seurat(seurat_integrate)

ider <- getIDEr(seu.integrated, batch.by.var = "orig.ident", downsampling.size = 11, verbose = TRUE)
print("finished!")
seu.integrated <- estimateProb(seu.integrated, ider)
p1 <- scatterPlot(seu.integrated, "tsne", colour.by = "similarity")
p2 <- scatterPlot(seu.integrated, "tsne", colour.by = "pvalue") 
plot_grid(p1,p2, ncol = 2)