library(ggplot2)
library(grid)
library(gridExtra)
library(glue)
library(dplyr)
library(Seurat)
library(sva)
library(stringr)
library(ggpca)
library(Hotelling)

# supercells <- readRDS("../SuperCell/supercells/5_10_supercell.rds")
# 
# # Microarray
# microarray <- read.delim('osd_324/GLDS-324_microarray_2015-10-07_Pajevic_RMA-GENE-FULL_Group.txt', )
# metadata <- read.delim2('osd_324/s_GLDS-324.txt')
# 
# colnames(metadata)[2] <- "sample"
# colnames(metadata)[16] <- "gravity"
# colnames(metadata)[19] <- "time"
# 
# metadata$sample <- sub(" ", "_", str_c(metadata$gravity, metadata$time, sep = "_"))
# 
# microarray$Gene.Symbol <- trimws(microarray$Gene.Symbol)
# colnames(microarray)[2:13] <- metadata$sample
# 
# microarray <- microarray[!microarray$Gene.Symbol == "---",]
# microarray <- microarray[!is.na(microarray$Gene.Symbol),]
# microarray <- microarray[!duplicated(microarray$Gene.Symbol),]
# 
# rownames(microarray) <- microarray$Gene.Symbol
# 
# micro_meta <- metadata
# micro_meta <- micro_meta[,c(2,16,19)]
# micro_meta$ligand <- rep("NO", 12)
# micro_meta$type <- rep("Microarray", 12)
# 
# microarray <- microarray[2:13]
# microarray_mat <- as.matrix(microarray)
# 
# microarray <- CreateSeuratObject(counts=microarray, assay = "microarray", meta.data = micro_meta)
# microarray <- NormalizeData(microarray)
# microarray <- FindVariableFeatures(microarray)
# microarray <- ScaleData(microarray)
# microarray <- RunPCA(microarray, npcs = 11)
# microarray <- RunUMAP(microarray, dims = 1:11, n.neighbors = 5, reduction = "pca", return.model = TRUE)
# saveRDS(microarray, "microarray_intersect.rds")
# 
# 
# # Proteomics
# proteomics <- readRDS("day1_3_7.rds")
# time <- integer(length(colnames(proteomics)))
# ligand <- character(length(colnames(proteomics)))
# gravity <- character(length(colnames(proteomics)))
# sample <- colnames(proteomics)
# 
# time[grep("Day.1", colnames(proteomics))] <- 1
# time[grep("Day.3", colnames(proteomics))] <- 3
# time[grep("Day.7", colnames(proteomics))] <- 7
# 
# ligand[grep("No", colnames(proteomics))] <- "NO"
# ligand[grep("with", colnames(proteomics))] <- "WITH"
# 
# gravity[grep("Micro", colnames(proteomics))] <- "MICRO"
# gravity[grep("Normal", colnames(proteomics))] <- "NORMAL"
# 
# sample <- as.factor(gsub(".{2}$", "", sample))
# 
# proteomics_meta <- data.frame(time, ligand, gravity, sample)
# proteomics_meta$type <- rep("Proteomics", 47)
# 
# rownames(proteomics_meta) <- colnames(proteomics)
# proteomics_multi <- proteomics
# 
# proteomics <- CreateSeuratObject(counts=proteomics, assay = "proteomics", meta.data = proteomics_meta)
# 
# # Integration using comBat Seq
# 
# proteomics_multi <- as.data.frame(proteomics_multi)
# microarray_mat <- as.data.frame(microarray_mat)
# 
# common_genes <- intersect(rownames(proteomics_multi), rownames(microarray_mat))
# 
# # Subset both data frames to only those genes
# proteomics_common <- proteomics_multi[common_genes, , drop = FALSE]
# microarray_common <- microarray_mat[common_genes, , drop = FALSE]
# 
# # Combine them side-by-side
# merged_data <- cbind(proteomics_common, microarray_common)
# merged_meta <- rbind(proteomics_meta, micro_meta)
# #
# merged_meta$time <- as.factor(merged_meta$time)
# merged_meta$gravity <- as.factor(merged_meta$gravity)
# merged_meta$sample <- as.factor(merged_meta$sample)
# merged_meta$type <- as.factor(merged_meta$type)
# merged_meta$ligand <- as.factor(merged_meta$ligand)
# 
# combat2 <- ComBat(merged_data, merged_meta$type)
# 
# combat <- CreateSeuratObject(counts=combat2, assay = "combat", meta.data = merged_meta)
# 
# combat@meta.data$orig.ident <- rep("ComBat", 59)
# combat@meta.data$time <- (merged_meta$time)
# combat@meta.data$gravity <- (merged_meta$gravity)
# combat@meta.data$sample <- (merged_meta$sample)
# combat@meta.data$type <- (merged_meta$type)
# combat@meta.data$ligand <- (merged_meta$ligand)


integrate_and_transfer <- function(supercells, proteomics, filename, pcs, n.neighbors, k.anchor, k.weight, k.score) {

  intersection <- intersect(rownames(supercells$seuratMC), rownames(proteomics))
  print(length(intersection))
  
  proteomics <- subset(proteomics, features = intersection)
  proteomics <- NormalizeData(proteomics)
  proteomics <- FindVariableFeatures(proteomics)
  proteomics <- ScaleData(proteomics)
  proteomics <- RunPCA(proteomics, npcs = pcs)
  proteomics <- RunUMAP(proteomics, dims = 1:pcs, n.neighbors = n.neighbors, reduction = "pca", return.model = TRUE)
  saveRDS(proteomics, "proteomics_intersect.rds")
  
  supercells <- subset(supercells$seuratMC, features = intersection)
  
  # Seurat data transfer
  # transfer.anchors <- FindTransferAnchors(proteomics, supercells, reduction = "pcaproject",
  #                                         reference.reduction = "pca", k.score = k.score,
  #                                         project.query = FALSE, dims = 1:pcs, verbose = TRUE)
  # 
  # transfer.data <- TransferData(anchorset = transfer.anchors, refdata = proteomics$sample, reference = proteomics,
  #                               query = supercells, k.weight = k.weight, weight.reduction = "pcaproject")
  # 
  # integrate.embeddings <- IntegrateEmbeddings(anchorset = transfer.anchors, reference = proteomics, query = supercells,
  #                                             new.reduction.name = "ref.pca", k.weight = k.weight)
  # 
  # map.query <- ProjectUMAP(query = integrate.embeddings, query.reduction= "ref.pca", reference = proteomics,
  #                          reference.reduction = "pca", reduction.model = "umap")
  # 
  # testplot <- ProjectDimReduc(supercells, proteomics, mode = "pcaproject", reference.reduction = "pca", reduction.name = "pca_proj")
  # 
  # pdf(file = glue("seurat_transfer_plots/{filename}_pca.pdf"), width=16, height=9)
  # prot <- DimPlot(proteomics, reduction = "pca", group.by = "ligand") + xlim(-40, 40) + ylim(-30, 20)
  # sc <- DimPlot(map.query, reduction = "ref.pca", group.by = "time")+ xlim(-40, 40) + ylim(-30, 20)
  # grid.arrange(prot, sc, nrow=1, top="pca (prot + sc) - projectUmap")
  # dev.off()
  # 
  # pdf(file = glue("seurat_transfer_plots/{filename}_umap.pdf"), width=16, height=9)
  # prot <- DimPlot(proteomics, reduction = "umap", group.by = "ligand") + xlim(-10, 10) + ylim (-8, 8)
  # sc <- DimPlot(map.query, reduction = "ref.umap", group.by = "time")+ xlim(-10, 10) + ylim (-8, 8)
  # grid.arrange(prot, sc, nrow=1, top="umap (prot + sc)")
  # dev.off()
  # 
  # pdf(file = glue("seurat_transfer_plots/{filename}_pcatest.pdf"), width=16, height=9)
  # prot <- DimPlot(proteomics, reduction = "pca", group.by = "ligand") + xlim(-40, 50) + ylim(-30, 20)
  # sc <- DimPlot(testplot, reduction = "pca_proj", group.by = "time") + xlim(-40, 50) + ylim(-30, 20)
  # grid.arrange(prot, sc, nrow=1, top="pca (prot + sc) - projectDimReduc")
  # dev.off()
  
  # Seurat data integration
  seurat.list <- list(proteomics, supercells)
  
  features <- SelectIntegrationFeatures(object.list = seurat.list)
  
  print("features")
  
  integration.anchors <- FindIntegrationAnchors(object.list = seurat.list, 
                                                anchor.features = features, dims = 1:pcs,
                                                reduction = "cca", k.anchor = k.anchor)
  print("integration.anchors")
  
  seurat.combined <- IntegrateData(anchorset = integration.anchors, k.weight = k.weight)
  
  DefaultAssay(seurat.combined) <- "integrated"
  
  seurat.combined <- ScaleData(seurat.combined, verbose = FALSE)
  seurat.combined <- RunPCA(seurat.combined, npcs = pcs, verbose = FALSE)
  seurat.combined <- RunUMAP(seurat.combined, reduction = "pca", dims = 1:pcs)
  seurat.combined <- FindNeighbors(seurat.combined, reduction = "pca", dims = 1:pcs)
  seurat.combined <- FindClusters(seurat.combined, resolution = 0.5)

  pdf(file = glue("seurat_integrate_plots/{filename}_pca.pdf"), width=16, height=9)
  sc <- DimPlot(seurat.combined, reduction = "pca", group.by = "sample", split.by = "orig.ident")
  grid.arrange(sc)
  dev.off()
  
  pdf(file = glue("seurat_integrate_plots/{filename}_ligand_split_pca.pdf"), width=16, height=9)
  sc <- DimPlot(seurat.combined, reduction = "pca", group.by = "sample", split.by = "ligand")
  grid.arrange(sc)
  dev.off()
  
  pdf(file = glue("seurat_integrate_plots/{filename}_umap.pdf"), width=16, height=9)
  sc <- DimPlot(seurat.combined, reduction = "umap", group.by = "sample", split.by = "orig.ident")
  grid.arrange(sc)
  dev.off()
  
  pdf(file = glue("seurat_integrate_plots/{filename}_sample.pdf"), width=16, height=9)
  sc <- DimPlot(seurat.combined, reduction = "umap", group.by = "sample", split.by = "orig.ident")
  grid.arrange(sc)
  dev.off()
  
  return(c(proteomics, supercells, seurat.combined))
}

# combatobj <- integrate_and_transfer(supercells = supercells, proteomics = combat,
#                                     n.neighbors = 24, k.anchor = 27,
#                                     pcs = 20, k.weight = 27, k.score = 27, filename = "combat")
# testobj <- integrate_and_transfer(supercells = supercells, proteomics = proteomics,
#                                   n.neighbors = 10, k.anchor = 27,
#                                   pcs = 20, k.weight = 27, k.score = 27, filename = "proteomics")
# saveRDS(testobj[[3]], file = "seurat.combined.rds")
# #

# 
# saveRDS(combatobj[[3]], file = "combat.combined.rds")


# pdf(file = glue("dim_reduct_plots/ComBat_PCA.pdf"), width=16, height=9)
# prot <- DimPlot(combatobj[[1]], reduction = "pca", group.by = "sample", label = TRUE)
# grid.arrange(prot, nrow=1, top="ComBat PCA")
# dev.off()
# 
# pdf(file = glue("dim_reduct_plots/Integrated_UMAP.pdf"), width=16, height=9)
# prot <- DimPlot(combatobj[[3]], reduction = "umap", group.by = "sample", split.by = "orig.ident")
# grid.arrange(prot, nrow=1, top="ComBat PCA")
# dev.off()

pca_func <- function(proteomics_intersect, group.by, reduction = "pca",
                     split.by = NULL,
                     colours = NULL,
                     legend_title = NULL,
                     match_border = FALSE,
                     ellipse_expand = 1.1  # smaller = skinnier ellipse
) {
  library(rlang)
  library(ggplot2)
  library(patchwork)
  library(scales)
  library(dplyr)
  
  if (!is.character(group.by) || length(group.by) != 1) {
    stop("`group.by` must be a single column name string, e.g. 'plot_pca_1'")
  }
  
  pca_var <- proteomics_intersect[["pca"]]@stdev
  pca_percent <- round(100 * (pca_var^2 / sum(pca_var^2)), 1)
  
  plots <- DimPlot(proteomics_intersect,
                   reduction = reduction,
                   group.by = group.by,
                   split.by = split.by,
                   combine = FALSE)
  
  plots <- lapply(plots, function(p) {
    df <- p$data
    grp_levels <- levels(factor(df[[group.by]]))
    
    # colour setup
    if (is.null(colours)) {
      pal <- hue_pal()(length(grp_levels))
      colours_vec <- setNames(pal, grp_levels)
    } else {
      if (is.null(names(colours))) {
        if (length(colours) < length(grp_levels)) stop("Not enough colours provided.")
        colours_vec <- setNames(colours[seq_along(grp_levels)], grp_levels)
      } else {
        missing_levels <- setdiff(grp_levels, names(colours))
        if (length(missing_levels) > 0)
          stop("Provided `colours` is missing levels: ", paste(missing_levels, collapse = ", "))
        colours_vec <- colours[grp_levels]
      }
    }
    
    # normal ellipses for e3 points
    if (match_border) {
      p <- p +
        stat_ellipse(
          data = df |> group_by(!!sym(group.by)) |> filter(n() >= 3),
          aes(x = PC_1, y = PC_2, fill = !!sym(group.by), colour = !!sym(group.by)),
          geom = "polygon", alpha = 0.25, level = 0.9, show.legend = FALSE
        )
    } else {
      p <- p +
        stat_ellipse(
          data = df |> group_by(!!sym(group.by)) |> filter(n() >= 3),
          aes(x = PC_1, y = PC_2, fill = !!sym(group.by)),
          colour = "black", geom = "polygon", alpha = 0.25, level = 0.9, show.legend = FALSE
        )
    }
    
    # --- draw ellipse using two points as foci ---
    two_point_groups <- df |>
      group_by(!!sym(group.by)) |>
      filter(n() == 2) |>
      pull(!!sym(group.by)) |>
      unique()
    
    for (grp in two_point_groups) {
      sub_df <- df[df[[group.by]] == grp, c("PC_1", "PC_2")]
      
      # foci
      f1 <- as.numeric(sub_df[1, ])
      f2 <- as.numeric(sub_df[2, ])
      c_dist <- sqrt(sum((f1 - f2)^2)) / 2  # c
      center <- (f1 + f2) / 2               # midpoint
      
      a <- ellipse_expand * c_dist
      b <- sqrt(a^2 - c_dist^2)
      
      # angle of major axis (radians)
      theta <- atan2(f2[2] - f1[2], f2[1] - f1[1])
      
      # parametric ellipse centered at origin
      t <- seq(0, 2 * pi, length.out = 300)
      ellipse_base <- cbind(a * cos(t), b * sin(t))
      
      # rotate and translate correctly
      R <- matrix(c(cos(theta), -sin(theta),
                    sin(theta),  cos(theta)), ncol = 2)
      ellipse_rot <- ellipse_base %*% R
      ellipse_df <- data.frame(
        PC_1 = ellipse_rot[, 1] + center[1],
        PC_2 = ellipse_rot[, 2] + center[2],
        group = grp
      )
      
      fill_col <- colours_vec[grp]
      border_col <- if (match_border) fill_col else "black"
      
      p <- p +
        geom_polygon(
          data = ellipse_df,
          aes(x = PC_1, y = PC_2),
          fill = fill_col,
          colour = border_col,
          alpha = 0.25,
          inherit.aes = FALSE,
          show.legend = FALSE
        )
    }
    
    p <- p +
      scale_colour_manual(name = legend_title %||% group.by,
                          values = colours_vec) +
      scale_fill_manual(values = colours_vec, guide = "none") +
      labs(
        x = paste0("PC_1 (", pca_percent[1], "% variance)"),
        y = paste0("PC_2 (", pca_percent[2], "% variance)")
      ) +
      theme_classic() 
    
    return(p)
  })
  
  wrap_plots(plots)
}

# plotprot <- str_c(testobj[[1]]$ligand, "ligand", "Day", testobj[[1]]$time, sep = "_")
# 
# testobj[[1]]@meta.data$plotprot <- plotprot
# 
# test2 <- pca_func(microarray, "sample", match_border = TRUE, ellipse_expand = 1.1)
# test3 <- pca_func(testobj[[1]], "plotprot", match_border = TRUE)
# test4 <- pca_func(testobj[[1]], "ligand", match_border = TRUE)
# test5 <- pca_func(testobj[[1]], "time", match_border = TRUE)
# test6 <- pca_func(testobj[[1]], "gravity", match_border = TRUE)

# combined_plot <- patchwork::wrap_plots(umap_plot)
# 
# testobj[[3]]$plot_umap_1 <- testobj[[3]]@meta.data %>% mutate(
#   composite = case_when(
#     orig.ident == "SeuratProject" ~ time,
#     TRUE ~ time
#     #TRUE ~ paste0(time, "_", ligand)
#   )) %>% pull(composite)
# 
# umap_plot <- DimPlot(testobj[[3]], reduction="umap", split.by = "orig.ident", 
#                      group.by = "plot_umap_1")
# umap_plot <- lapply(umap_plot, function(p) { p + stat_ellipse(aes(x = umap_1, y = umap_2, color = plot_umap_1), level = 0.9) }) 
# 
# combined_plot_2 <- patchwork::wrap_plots(umap_plot)


combatobj[[3]]$plot_umap_1 <- combatobj[[3]]@meta.data %>% mutate(
  composite = case_when(
    orig.ident != "ComBat" ~ time,
    type == "Proteomics" ~ ligand,
    type == "Microarray" ~ gravity
  )) %>% pull(composite)

test7 <- pca_func(combatobj[[3]], group.by = "gravity", split.by="orig.ident", match_border = TRUE)
test8 <- pca_func(combatobj[[3]], group.by = "gravity", split.by="orig.ident", 
                  reduction = "umap", match_border = TRUE)

pca_plot <- DimPlot(combatobj[[3]], reduction="pca", split.by = "orig.ident", 
                     group.by = "plot_umap_1")
pca_plot <- lapply(pca_plot, function(p) { p + stat_ellipse(aes(x = PC_1, y = PC_2, color = plot_umap_1), level = 0.9) }) 

umap_plot <- DimPlot(combatobj[[3]], reduction="umap", split.by = "orig.ident", 
                     group.by = "plot_umap_1")
umap_plot <- lapply(umap_plot, function(p) { p + stat_ellipse(aes(x = umap_1, y = umap_2, color = plot_umap_1), level = 0.9) }) 



# library(Hotelling)
# 
# k <- 20
# pca_df <- as.data.frame(Embeddings(testobj[[1]], "pca"))[, 1:k]
# pca_df$ligand <- testobj[[1]]@meta.data$ligand
# 
# pca_df %>% group_by(ligand) %>%
#   mutate(Dist = colMeans(as.matrix(dist(cbind(PC_1, PC_2)))))
# 
# x1 <- pca_df[pca_df$ligand == "NO", 1:k]
# x2 <- pca_df[pca_df$ligand == "WITH", 1:k]
# compute Hotelling's T2 -> F as in your previous conversion
hotelling_centroid_test <- function(x1, x2) {
  n1 <- nrow(x1); n2 <- nrow(x2); p <- ncol(x1)
  if (n1 < 2 || n2 < 2) stop("Each group needs >= 2 observations.")
  m1 <- colMeans(x1); m2 <- colMeans(x2)
  S1 <- cov(x1); S2 <- cov(x2)
  Spooled <- ((n1 - 1)*S1 + (n2 - 1)*S2) / (n1 + n2 - 2)
  diff <- m1 - m2
  T2 <- (n1 * n2) / (n1 + n2) * as.numeric(t(diff) %*% solve(Spooled) %*% diff)
  df1 <- p
  df2 <- n1 + n2 - p - 1
  if (df2 <= 0) stop("n1 + n2 - p - 1 must be > 0 for F conversion. Reduce p.")
  Fval <- ((n1 + n2 - p - 1) * T2) / (p * (n1 + n2 - 2))
  # get log p (upper tail). This will usually be finite and informative.
  logp <- pf(Fval, df1, df2, lower.tail = FALSE, log.p = TRUE)
  list(T2 = T2, F = Fval, df = c(df1, df2), logp = logp,
       p_estimate = if (is.finite(logp)) exp(logp) else 0)
}

# res <- hotelling_centroid_test(x1, x2)
# cat("T2 =", res$T2, "F =", res$F, "df =", paste(res$df, collapse = ","), "\n")
# cat("log10 p (approx) =", as.numeric(res$logp / log(10)), "\n")
# # if res$logp is -Inf, fallback to permutations (below)
# 
# res1 <- hotelling.test(x1, x2)
# p_one_sided <- res1$pval /2  # approximate one-sided p
# p_adj <- p.adjust(p_one_sided, method = "bonferroni")
# print(p_adj)
# 
# pca_df <- as.data.frame(Embeddings(testobj[[1]], "pca"))[, 1:k]
# pca_df$gravity <- testobj[[1]]@meta.data$gravity
# 
# x1 <- pca_df[pca_df$gravity == "MICRO", 1:k]
# x2 <- pca_df[pca_df$gravity == "NORMAL", 1:k]
# 
# res <- hotelling.test(x1, x2)
# p_one_sided <- res$pval /2  # approximate one-sided p
# p_adj <- p.adjust(p_one_sided, method = "bonferroni")
# print(p_adj)
# 
# pca_df <- as.data.frame(Embeddings(testobj[[1]], "pca"))[, 1:k]
# pca_df$time <- testobj[[1]]@meta.data$time
# 
# x1 <- pca_df[pca_df$time == "1", 1:k]
# x3 <- pca_df[pca_df$time == "3", 1:k]
# x7 <- pca_df[pca_df$time == "7", 1:k]
# 
# res <- hotelling.test(x1, x3)
# p_one_sided <- res$pval /2  # approximate one-sided p
# p_adj <- p.adjust(p_one_sided, method = "bonferroni")
# print(p_adj)
# 
# res <- hotelling.test(x7, x3)
# p_one_sided <- res$pval /2  # approximate one-sided p
# p_adj <- p.adjust(p_one_sided, method = "bonferroni")
# print(p_adj)
# 
# res <- hotelling.test(x1, x7)
# p_one_sided <- res$pval /2  # approximate one-sided p
# p_adj <- p.adjust(p_one_sided, method = "bonferroni")
# print(p_adj)
# 
# prot_list <- list(test4, test6, test5)
# pdf(file = "thesis_plots/prot_pca_2.pdf", height = 4, width = 13)
# grid.arrange(grobs = prot_list, nrow=1)
# dev.off()
# 
# k <- 20
# levels <- unique(testobj[[1]]$plotprot)
# pca_df <- as.data.frame(Embeddings(testobj[[1]], "pca"))[, 1:k]
# pca_df$plotprot <- testobj[[1]]@meta.data$plotprot
# 
# 
# 
# for (x in sort(levels)) {
#   test <- pca_df[pca_df$plotprot == x, 1:20]
#   
#   print(mean(dist(test)))
# }
# 
# distance <- c(34.99688, 28.98478, 40.74593, 12.94775, 12.3111,16.47107)
# 
# tabley <- data.frame(sort(levels), distance)
# grid.arrange(tableGrob(tabley))
# k <- 20
# combatobj[[1]]$plot_1 <- combatobj[[1]]@meta.data %>% mutate(
#   composite = case_when(
#     type == "Microarray" ~ "Cluster 1",
#     ligand == "WITH" ~ "Cluster 1",
#     TRUE ~ "Cluster 2"
#   )) %>% pull(composite)
# 
# test8 <- pca_func(combatobj[[1]], "plot_1", match_border = TRUE)
# 
# levels <- unique(combatobj[[1]]$plot_1)
# pca_df <- as.data.frame(Embeddings(combatobj[[1]], "pca"))[, 1:k]
# pca_df$plotprot <- combatobj[[1]]@meta.data$plot_1
# 
# for (x in sort(levels)) {
#   test <- pca_df[pca_df$plotprot == x, 1:20]
# 
#   print(mean(dist(test)))
# }
