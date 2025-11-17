library(Seurat)
library(glue)
library(ggplot2)
#library(ggpubr)
#library(patchwork)
#library(VennDiagram)
library(dplyr)
library(SuperCell)
library(grid)

# singlecells <- readRDS("Ost_all.rds")
# supercells <- readRDS("supercells.rds")

# plot2 <- ggplot(supercells$seuratMC@meta.data, aes(x=as.factor(time), y=purity)) + geom_boxplot(fill="slateblue", alpha=0.2) + xlab("sample")
# png(file=glue("plot/supercell_purity_time.png"), width=1080, height=1080)
# grid.draw(plot2)
# dev.off()
     
plot3 <- ggplot(supercells$seuratMC@meta.data, aes(x=as.factor(sample), y=size)) + geom_boxplot(outlier.shape=NA, fill="slateblue", alpha=0.2) + xlab("sample") + ylim(0, 100)
png(file=glue("plot/supercell_size_sample.png"), width=1080, height=1080)
grid.draw(plot3)
dev.off()
     
plot3 <- ggplot(supercells$seuratMC@meta.data, aes(x=as.factor(sample), y=separation)) + geom_boxplot(outlier.shape=1, fill="slateblue", alpha=0.2) + xlab("sample")
png(file=glue("plot/supercell_separation_sample.png"), width=1080, height=1080)
grid.draw(plot3)
dev.off()
     
plot3 <- ggplot(supercells$seuratMC@meta.data, aes(x=as.factor(sample), y=compactness)) + geom_boxplot(outlier.shape=1, fill="slateblue", alpha=0.2) + xlab("sample")
png(file=glue("plot/supercell_compactness_sample.png"), width=1080, height=1080)
grid.draw(plot3)
dev.off()

plot3 <- ggplot(supercells$seuratMC@meta.data, aes(x=as.factor(sample), y=INV)) + geom_boxplot(outlier.shape=1, fill="slateblue", alpha=0.2) + xlab("sample") 
png(file=glue("plot/supercell_INV_sample.png"), width=1080, height=1080)
grid.draw(plot3)
dev.off()
     
