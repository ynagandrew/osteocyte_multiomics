library(stats)
library(limma)
library(grid)
library(glue)
library(forcats)
library(readxl)
library(clusterProfiler)
library(org.Mm.eg.db)
library(msigdbr)
library(enrichplot)
library(stringr)
library(gridExtra)
library(ggplot2)

osd_324 <- read_xlsx("OSD-324/GLDS-324_microarray_2015-10-07_Pajevic_analysis.xlsx")
osd_635 <- read.csv("OSD-635/GLDS-608_rna_seq_differential_expression_GLbulkRNAseq.csv")
agoro <- read.csv("agoro/data_sheet_1.csv")
kaur <- read_xlsx("kaur/can-24-0857_supplementary_table_1_suppst1.xlsx")
wang <- read_xlsx("wang/41467_2021_26571_MOESM9_ESM.xlsx")

youlten <- read.csv("youlten/Gene_annotation_w_activity_enrichment_signature.csv")

# ost_trans_sig <- filter(youlten, Osteocyte_transcriptome_signature == TRUE)
# 
# ost_trans_sig$EntrezID
# 
# go_test <- enrichGO(gene = ost_trans_sig$EntrezID,
#                universe = youlten$EntrezID,
#                OrgDb = org.Mm.eg.db,
#                pvalueCutoff = 0.05,
#                pAdjustMethod = "BH",
#                keyType = "ENTREZID")