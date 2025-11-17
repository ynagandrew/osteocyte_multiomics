library(dplyr)
library(limma)
library(clusterProfiler)
library(enrichplot)
library(org.Mm.eg.db)
library(glue)
library(readxl)
library(gridExtra)
library(GOSemSim)
library(UCell)

library(ggplot2)
library(ggpubr)

organism = org.Mm.eg.db

bpGO <- godata(organism, ont="BP")
ccGO <- godata(organism, ont="CC")

ost_micro_sig <- readRDS("../osteocyte_microgravity_signature_name.rds")

normd <- read.delim('OSD-324/GLDS-324_microarray_2015-10-07_Pajevic_RMA-GENE-FULL_Group.txt',
                     )
metadata <- read.delim2('OSD-324/s_GLDS-324.txt')

normd$Gene.Symbol <- trimws(normd$Gene.Symbol)
colnames(normd)[2:13] <- metadata$Sample.Name

gravity <- as.factor(metadata$Factor.Value.Spaceflight)
time <- as.factor(metadata$Factor.Value.Time)

csv <- normd[2:13]
design <- model.matrix(~gravity + time)

individual_comparison <- function(data, design, contrasts) {

  #TOPLOT: Variability vs Count Size
  v <- voom(data, design, plot=FALSE)

  fit <- lmFit(v, design)
  coefs <- (fit$coefficients)
  fit <- contrasts.fit(fit, contrasts)
  fit <- eBayes(fit)

  return (fit)
}

get_coefs <- function(data, design) {

  #TOPLOT: Variability vs Count Size
  v <- voom(data, design, plot=FALSE)

  fit <- lmFit(v, design)
  fit <- eBayes(fit)

  return (fit)
}


# DEGenes Data

osd_324 <- read_xlsx("OSD-324/GLDS-324_microarray_2015-10-07_Pajevic_analysis.xlsx")

modify_csv <- function(csv) {
  colnames(csv) <- csv[1,]
  csv <- csv[-1,]

  csv <- csv %>% filter(`FDR q (filtered)` > 0.05)

  gene_list <- csv$`fold change`
  gene_list <- as.numeric(gene_list)
  names(gene_list) <- csv$Symbol

  gene_list <- gene_list[!names(gene_list) %in% c("---")]
  gene_list <- na.omit(gene_list)
  gene_list <- gene_list[!is.na(names(gene_list))]
  gene_list <- gene_list[!duplicated(names(gene_list))]

  gene_list <- sort(gene_list, decreasing = TRUE)

  return(gene_list)
}

gse_plots <- function(gse, filename) {
  pdf(file=glue("gse_plots/{filename}.pdf"), width=16, height=9)
  grid.arrange(dotplot(gse, showCategory = 5, split=".sign"), nrow=1, top=glue("{filename} GSE"))
  dev.off()
}

four_v_two <- osd_324[c(2,32:36)]
four_v_two <- modify_csv(four_v_two)

six_v_two <- osd_324[c(2,37:41)]
six_v_two <- modify_csv(six_v_two)

six_v_four <- osd_324[c(2,42:46)]
six_v_four <- modify_csv(six_v_four)

day2 <- osd_324[c(2,47:51)]
day2 <- modify_csv(day2)
day4 <- osd_324[c(2,52:56)]
day4 <- modify_csv(day4)
day6 <- osd_324[c(2,57:61)]
day6 <- modify_csv(day6)

get_genes <- function(listy) {
  outs_genes = c()
  for (csv in listy) {
    genes = names(csv)[1:50]
    genes_end <- tail(names(csv), n=50)

    outs_genes <- union(outs_genes, genes)
    outs_genes <- union(outs_genes, genes_end)
    print(length(outs_genes))
  }
  return(outs_genes)
}

genes <- list(four_v_two, six_v_two, six_v_four, day2, day4, day6)

all_genes <- get_genes(genes)

test <- intersect(all_genes, ost_micro_sig)
print(test)

gse_2_CC <- gseGO(geneList=day2,
                 ont ="CC",
                 keyType = "SYMBOL",
                 minGSSize = 3,
                 maxGSSize = 800,
                 pvalueCutoff = 0.05,
                 verbose = TRUE,
                 OrgDb = organism,
                 pAdjustMethod = "BH")

gse_4_CC <- gseGO(geneList=day4,
               ont ="CC",
               keyType = "SYMBOL",
               minGSSize = 3,
               maxGSSize = 800,
               pvalueCutoff = 0.05,
               verbose = TRUE,
               OrgDb = organism,
               pAdjustMethod = "BH")

gse_6_CC <- gseGO(geneList=day6,
                  ont ="CC",
                  keyType = "SYMBOL",
                  minGSSize = 3,
                  maxGSSize = 800,
                  pvalueCutoff = 0.05,
                  verbose = TRUE,
                  OrgDb = organism,
                  pAdjustMethod = "BH")

gse_2_BP <- gseGO(geneList=day2,
               ont ="BP",
               keyType = "SYMBOL",
               minGSSize = 3,
               maxGSSize = 800,
               pvalueCutoff = 0.05,
               verbose = TRUE,
               OrgDb = organism,
               pAdjustMethod = "BH")

gse_4_BP <- gseGO(geneList=day4,
               ont ="BP",
               keyType = "SYMBOL",
               minGSSize = 3,
               maxGSSize = 800,
               pvalueCutoff = 0.05,
               verbose = TRUE,
               OrgDb = organism,
               pAdjustMethod = "BH")

gse_6_BP <- gseGO(geneList=day6,
                  ont ="BP",
                  keyType = "SYMBOL",
                  minGSSize = 3,
                  maxGSSize = 800,
                  pvalueCutoff = 0.05,
                  verbose = TRUE,
                  OrgDb = organism,
                  pAdjustMethod = "BH")

gse_plots(gse_2_CC, "2_gravs_CC")
gse_plots(gse_4_CC, "4_gravs_CC")
gse_plots(gse_6_CC, "6_gravs_CC")

gse_plots(gse_2_BP, "2_gravs_BP")
gse_plots(gse_4_BP, "4_gravs_BP")
gse_plots(gse_6_BP, "6_gravs_BP")

BP_2 <- gse_2_BP@result
BP_2 <- BP_2[with(BP_2, order(-BP_2$NES)), ]$ID[1:50]

BP_4 <- gse_4_BP@result
BP_4 <- BP_4[with(BP_4, order(-BP_4$NES)), ]$ID[1:50]

BP_6 <- gse_6_BP@result
BP_6 <- BP_6[with(BP_6, order(-BP_6$NES)), ]$ID[1:50]

CC_2 <- gse_2_CC@result
CC_2 <- CC_2[with(CC_2, order(-CC_2$NES)), ]$ID[1:50]

CC_4 <- gse_4_CC@result
CC_4 <- CC_4[with(CC_4, order(-CC_4$NES)), ]$ID[1:50]

CC_6 <- gse_6_CC@result
CC_6 <- CC_6[with(CC_6, order(-CC_6$NES)), ]$ID[1:50]

# micro_sig_CC <- cc@result[1:150,]$ID
# 
# micro_sig_BP <- bp@result[1:150,]$ID

print(mgoSim(BP_2, BP_4, bpGO))
# print(mgoSim(BP_2, BP_6, bpGO))
# print(mgoSim(BP_4, BP_6, bpGO))

print(mgoSim(CC_2, CC_4, ccGO))
# print(mgoSim(CC_2, CC_6, ccGO))
# print(mgoSim(CC_4, CC_6, ccGO))

# ontology_sig_CC <- union(CC_2, union(CC_4, CC_6))
# ontology_sig_BP <- union(BP_2, union(BP_4, BP_6))
# 
# print("here we go!")
# print(mgoSim(micro_sig_CC, ontology_sig_CC, ccGO))
# print(mgoSim(micro_sig_BP, ontology_sig_BP, bpGO))

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

ost_micro_sig <- list(Ost_Micro = ost_micro_sig)

ucell <- ScoreSignatures_UCell(microarray, ost_micro_sig, maxRank = 5000)

metadata <- cbind(metadata, ucell)

png(file=glue("../thesis_plots/microarray_ucell.png"), width=720, height=1080)
plot3 <- ggplot(metadata, aes(x=as.factor(Factor.Value.Spaceflight.), y=Ost_Micro_UCell)) +
  geom_boxplot(outlier.shape = NA) + xlab("Microarray Gravity") + ylab("Osteocyte Microgravity Signature Score") +
  #stat_compare_means(size=10)
  geom_bracket(xmin='Ground Control', xmax='Space Flight', y.position=0.27, label = "Wilcoxon, p = 0.041, maxRank = 5000", label.size = 10)

plot3 <- plot3 + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                                    panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), text = element_text(size=30),
                                    axis.text = element_text(size=30))
grid.arrange(plot3, nrow=1)
dev.off()
