library(Seurat)
library(UCell)
library(gridExtra)
library(glue)
library(ggpubr)
source("/scratch/user/s4704475/SuperCell/mc_projection.R")

# symnum.args <- list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, Inf), symbols = c("****", "***", "**", "*", "ns"))
# 
# ost_micro_sig <- readRDS("osteocyte_microgravity_signature_name.rds")
# 
# allost <- readRDS("../SuperCell/Osteocyte_v2025.rds")
# supercells <- readRDS("../SuperCell/supercells/5_10_supercell.rds")
# super_5_20 <- readRDS("../SuperCell/supercells/5_20_supercell.rds")
# 
# super_5_250 <- readRDS("../SuperCell/supercells/5_250_supercell.rds")
# super_100_250 <- readRDS("../SuperCell/supercells/100_250_supercell.rds")
#
# ost_micro_sig <- list(Ost_Micro = ost_micro_sig)
# 
# supercells$seuratMC <- AddModuleScore_UCell(obj = supercells$seuratMC, features = ost_micro_sig, maxRank = 2000)

png(file = "thesis_plots/metacell_a.png", width = 720, height = 720)
all <- DimPlot(allost, group.by = "time")
grid.arrange(all, nrow=1)
dev.off()

png(file = "thesis_plots/metacell_b.png", width = 720, height = 720)
all <- DimPlot(super_5_20$seuratMC, group.by = "time")
grid.arrange(all, nrow=1)
dev.off()

png(file = "thesis_plots/metacell_a.png", width = 720, height = 720)

grid.arrange(all, nrow=1)
dev.off()

png(file = "thesis_plots/metacell_b.png", width = 1920, height = 1080)
all1 <- DimPlot(allost, group.by = "time") + theme(text=element_text(size=28))
all <- DimPlot(super_5_20$seuratMC, group.by = "time")+ theme(text=element_text(size=28))
grid.arrange(all1, all, nrow=1)
dev.off()
# 
# png(file = glue("thesis_plots/metacell_b.png"), width=720, height=720)
# prot <- mc_projection(allost, supercells, metacell.label = "sample", sc.label="sample", pt_size=0.01)
# grid.arrange(prot, nrow=1)
# dev.off()
# 
# png(file = glue("thesis_plots/metacell_c.png"), width=720, height=720)
# prot <- mc_projection(allost, super_5_250, metacell.label = "sample", sc.label="sample", pt_size=0.01)
# grid.arrange(prot, nrow=1)
# dev.off()
# 
# png(file = glue("thesis_plots/metacell_d.png"), width=720, height=720)
# prot <- mc_projection(allost, super_100_250, metacell.label = "sample", sc.label="sample", pt_size=0.01)
# grid.arrange(prot, nrow=1)
# dev.off()
# 
# png(file = "thesis_plots/metacell_e.png", width = 720, height = 720)
# all <- FeaturePlot(supercells$seuratMC, features = "Ost_Micro_UCell")
# grid.arrange(all, nrow=1)
# dev.off()

# png(file=glue("thesis_plots/metacell_g.png"), width=1080, height=1080)
# plot3 <- ggplot(supercells$seuratMC@meta.data, aes(x=as.factor(sample), y=Ost_Micro_UCell)) + geom_boxplot(outlier.shape=NA, fill="slateblue", alpha=0.2) + xlab("Ost_Micro") + ylim(0.125, 0.375) + 
#   grid.arrange(plot3, nrow=1)
# dev.off()

png(file=glue("thesis_plots/metacell_f.png"), width=1080, height=1080)
plot3 <- ggplot(supercells$seuratMC@meta.data, aes(x=as.factor(stiffness), y=Ost_Micro_UCell)) +
  geom_boxplot(outlier.shape = NA) + xlab("Metacell Ligand Stiffness") + ylab("Osteocyte Microgravity Signature Score") +
  #stat_compare_means(size=10, symnum.args = symnum.args) +
  geom_bracket(xmin='L', xmax='S', y.position=0.275, label = "Wilcoxon, p < 2.2e-16", label.size = 10)

plot3 <- plot3 + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), text = element_text(size=30),
                       axis.text = element_text(size=30))
grid.arrange(plot3, nrow=1)
dev.off()


png(file=glue("thesis_plots/metacell_f.png"), width=1080, height=1080)
plot3 <- ggplot() + 
  geom_boxplot(aes(x="snRNASeq", y=allost$nFeature_RNA/21071)) +
  geom_boxplot(aes(x="Metacells", y=super_5_20$seuratMC$nFeature_RNA/21071)) +
  ylab("Density") + xlab("") + stat_compare_means() 
  #geom_bracket(xmin='L', xmax='S', y.position=0.275, label = "Wilcoxon, p < 2.2e-16", label.size = 10)

plot3 <- plot3 + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                                    panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), text = element_text(size=30),
                                    axis.text = element_text(size=30))
grid.arrange(plot3, nrow=1)
dev.off()

print(aggregate(supercells$seuratMC@meta.data$Ost_Micro_UCell, list(supercells$seuratMC@meta.data$ligand), mean))
