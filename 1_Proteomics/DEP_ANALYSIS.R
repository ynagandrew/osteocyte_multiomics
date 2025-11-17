# 1. Library Packages
library(DEP)
library(SummarizedExperiment)
library(dplyr)
library(rlang)

library(clusterProfiler)
library(org.Mm.eg.db)

library(glue)

library(enrichplot)
library(ggplot2)
library(ggpubr)
library(patchwork)
library(VennDiagram)
library(gridExtra)
library(ggplotify)
library(RColorBrewer)
library(viridis)

#library(wordcloud)

# 2. Load Datasets from CSV
day_1_no <- read.csv("day_1_no_ligand.csv")
day_1_with <- read.csv("day_1_with_ligand.csv")

day_3_no <- read.csv("day_3_no_ligand.csv")
day_3_with <- read.csv("day_3_with_ligand.csv")

day_7_no <- read.csv("day_7_no_ligand.csv")
day_7_with <- read.csv("day_7_with_ligand.csv")

# 2.1 Generate Venn Diagrams

# 3. Merge Datasets & QC

merge_dfs <- function(dfs) {
  if (length(dfs) == 2) {
    return(inner_join(dfs[[1]], dfs[[2]], by = "PG.UniProtIds"))
  } else {
    first_merge <- inner_join(dfs[[1]], dfs[[2]], by = "PG.UniProtIds")
    return(merge_dfs(c(list(first_merge), dfs[-c(1, 2)])))
  }
}

merge_dfs_full <- function(dfs) {
  if (length(dfs) == 2) {
    return(full_join(dfs[[1]], dfs[[2]], by = "PG.UniProtIds"))
  } else {
    first_merge <- full_join(dfs[[1]], dfs[[2]], by = "PG.UniProtIds")
    return(merge_dfs(c(list(first_merge), dfs[-c(1, 2)])))
  }
}

merge_proteomics_dfs <- function(df_merged) {
  
  df_merged[df_merged == "NaN"] <- NA
  
  # Handle day 7 dataset missing genes column
  if (is.null(df_merged$PG.Genes)) {
    data_unique <- make_unique(df_merged, "PG.Genes.y", "PG.UniProtIds", delim = ";")
  } else {
    data_unique <- make_unique(df_merged, "PG.Genes", "PG.UniProtIds", delim = ";")
  }
  
  LFQ_columns <- grep("_\\d+$", colnames(data_unique))
  data_se <- make_se_parse(data_unique, LFQ_columns)
  
  data_filt <- filter_missval(data_se, thr = 0)
  grid.draw(plot_detect(data_filt))
  
  data_norm <- normalize_vsn(data_filt)
  #grid.draw(plot_normalization(data_filt, data_norm))
  
  data_imp <- impute(data_norm, fun = "QRILC")
  #plot_imputation(data_norm, data_imp)
  data_diff_all_contrasts <- test_diff(data_imp, type = "all")
  dep <- add_rejections(data_diff_all_contrasts, alpha = 0.05, lfc = log2(1.5))
  
  return (dep)
}

# df1_full <- merge_dfs_full(list(day_1_no, day_1_with))
# df3_full <- merge_dfs_full(list(day_3_no, day_3_with))
# df3_full <- merge_dfs_full(list(day_3_no, day_3_with))
# df_all_full <- merge_dfs_full(list(df1_full, df3_full, df7_full))

df1 <- merge_dfs(list(day_1_no, day_1_with))
df3 <- merge_dfs(list(day_3_no, day_3_with))
df7 <- merge_dfs(list(day_7_no, day_7_with))
# df1_3_7 <- merge_dfs(list(df1, df3, df7))

# Initial analysis w/ DEP package (See limma.R for amended analysis)
# test1_full <- merge_proteomics_dfs(df1_full)
# test3_full <- merge_proteomics_dfs(df3_full)
# test7_full <- merge_proteomics_dfs(df7_full)
# test_all_full <- merge_proteomics_dfs(df_all_full)
# 
test1 <- merge_proteomics_dfs(df1)
test3 <- merge_proteomics_dfs(df3)
test7 <- merge_proteomics_dfs(df7)
#
# write.csv(df1, file = "day1.csv")
# write.csv(df3, file = "day3.csv")
# write.csv(df7, file = "day7.csv")
# 
# write.csv(df_all_full, file = "day1_3_7_full.csv")
# write.csv(df_1_3_7, file="day1_3_7.csv")


# 4. Generate proteomics plots

plot_proteomics <- function(df, day) {
  plot1 <- plot_pca(df, x = 1, y = 2, n = 500, point_size = 4)
  pdf(file=glue("plots/day_{day}_pca.pdf"))
  grid.draw(plot1)
  dev.off()
  
  conditions <- list(
    glue("Normal.G_Day.{day}_with.ligand__vs_Micro.G_Day.{day}_with.ligand_"),
    glue("Normal.G_Day.{day}_No.ligand__vs_Micro.G_Day.{day}_No.ligand_"),
    glue("Normal.G_Day.{day}_No.ligand__vs_Normal.G_Day.{day}_with.ligand_"),
    glue("Micro.G_Day.{day}_No.ligand__vs_Micro.G_Day.{day}_with.ligand_")
  )
  
  plots <- list()
  
  for (cond in conditions) {
    
    plots <- append(plots, list(as.grob(plot_volcano(df, contrast = cond, label_size = 2, add_names = FALSE))))
    
  }
  
  png(file=glue("plots/day_{day}_volcano.png"), width=1080, height=1080)
  combined_grobs <- grid.arrange(grobs = plots, ncol = 2)
  
  dev.off()
  return()
}


# 5. Generate Gene Set Enrichment Analysis & Plots

go_enrichment <- function(df, day) {
  df <- get_results(df)
    organism = org.Mm.eg.db
  
  # conditions <- list(
  #   glue("Normal.G_Day.{day}_with.ligand__vs_Micro.G_Day.{day}_with.ligand__ratio"),
  #   glue("Normal.G_Day.{day}_No.ligand__vs_Micro.G_Day.{day}_No.ligand__ratio")
  # )
  
  conditions <- list(
    glue("Normal.G_Day.{day}_No.ligand__vs_Normal.G_Day.{day}_with.ligand__ratio"),
    glue("Micro.G_Day.{day}_No.ligand__vs_Micro.G_Day.{day}_with.ligand__ratio")
  )
  
  dotplots <- list()
  networkplots <- list()
  heatplots <- list()
  
  gses <- character()
  
  for (cond in conditions) {
    df_input <- df[[cond]]
    names(df_input) <- df$ID
    gene_list<-na.omit(df_input)
    gene_list = sort(gene_list, decreasing = TRUE)
    gse <- gseGO(geneList=gene_list,
                 ont ="BP",
                 keyType = "UNIPROT",
                 minGSSize = 3,
                 maxGSSize = 800,
                 pvalueCutoff = 0.05,
                 verbose = TRUE,
                 OrgDb = organism,
                 pAdjustMethod = "none")
    
    require(DOSE)
    dotplots <- append(dotplots, list(as.grob(dotplot(gse, showCategory=5, split=".sign") + facet_grid(.~.sign) + labs(title = cond))))
    networkplots <- append(networkplots, list(as.grob(cnetplot(gse, foldChange = gene_list))))
    heatplots <- append(heatplots, list(as.grob(heatplot(gse, foldChange = gene_list, showCategory = 7))))
    
    gses <- c(gses, gse)
  }
  
  pdf(file=glue("thesis_plots/day_{day}_lig_go_dotplots.pdf"), width=9, height=11)
  combined_grobs <- grid.arrange(grobs = dotplots, ncol = 1, top=glue("Day {day} With Ligand v No Ligand"))
  
  dev.off()
  
  # png(file=glue("plots/day_{day}_go_networks.png"), width=1080, height=1920)
  # combined_grobs <- grid.arrange(grobs = networkplots, ncol = 2)
  # 
  # dev.off()
  # 
  # png(file=glue("plots/day_{day}_heatplot.png"), width=1920, height=1080)
  # combined_grobs <- grid.arrange(grobs = heatplots, ncol = 2)
  # 
  # dev.off()
  
  return (gses)
  
}

# 6. Generate KEGG Analysis & Plots

kegg_enrichment <- function(df, day) {
  df <- get_results(df)
  organism = "mmu"
  
  conditions <- list(
    glue("Normal.G_Day.{day}_with.ligand__vs_Micro.G_Day.{day}_with.ligand__ratio"),
    glue("Normal.G_Day.{day}_No.ligand__vs_Micro.G_Day.{day}_No.ligand__ratio"),
    glue("Normal.G_Day.{day}_No.ligand__vs_Normal.G_Day.{day}_with.ligand__ratio"),
    glue("Micro.G_Day.{day}_No.ligand__vs_Micro.G_Day.{day}_with.ligand__ratio")
  )
  
  for (cond in conditions) {
    df_input <- df[[cond]]
    names(df_input) <- df$ID
    gene_list<-na.omit(df_input)
    gene_list = sort(gene_list, decreasing = TRUE)
    kk2 <- gseKEGG(geneList     = gene_list,
                   organism     = "mmu",
                   minGSSize    = 3,
                   maxGSSize    = 800,
                   pvalueCutoff = 0.05,
                   pAdjustMethod = "BH",
                   keyType       = "uniprot")
    
    plot <- dotplot(kk2, showCategory=10, split=".sign") + labs(x = "enrichment distribution")
    pdf(file=glue("plots/{cond}_ridgeplot.pdf"))
    grid.draw(plot)
    dev.off()
  }
  
  return ()
  
}


# Functional Enrichment

gse1 <- go_enrichment(test1, "1")

gse3 <- go_enrichment(test3, "3")

gse7 <- go_enrichment(test7, "7")

gse_plots <- function(gses, day) {
  with <- dotplot(gses[[1]], showCategory=5, split=".sign") + labs(title=glue("Day {day}: Normal Gravity")) + facet_grid(.~.sign)
  no <- dotplot(gses[[2]], showCategory=5, split=".sign") + labs(title = "Microgravity") + facet_grid(.~.sign)
  pdf(file=glue("plots/day_{day}_go_dotplots.pdf"), width=9, height=13)
  combined_grobs <- grid.arrange(grobs = list(with, no), ncol = 1)
  dev.off()
}

gse_plots(gse1, "1")
gse_plots(gse3, "3")
gse_plots(gse7, "7")
