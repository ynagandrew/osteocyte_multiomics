library(Seurat)
library(philentropy)
library(ggplot2)
library(winch)
library(glue)
library(mixtools)
library(mclust)
library(dplyr)
library(gridExtra)
library(stringr)

#Using seurat.combined
seurat.combined <- readRDS("seurat.combined.rds")
DefaultAssay(seurat.combined) <- "integrated"

new.ident <- seurat.combined@meta.data$orig.ident

new.ident[new.ident == "Normal.G"] <- "Proteomics"
new.ident[new.ident == "Micro.G"] <- "Proteomics"

seurat.combined@meta.data$orig.ident <- new.ident
split <- SplitObject(seurat.combined, "orig.ident")

proteomics <- split[[1]]
seur_rna <- split[[2]]

# Using combat.combined
combat.combined <- readRDS("combat.combined.rds")
DefaultAssay(combat.combined) <- "integrated"

new.ident <- combat.combined@meta.data$orig.ident

combat.combined$plot_umap_1 <- combat.combined@meta.data %>% mutate(
  composite = case_when(
    orig.ident != "ComBat" ~ time,
    type == "Proteomics" ~ ligand,
    type == "Microarray" ~ gravity
  )) %>% pull(composite)

combat.combined@meta.data$orig.ident <- new.ident
split <- SplitObject(combat.combined, "orig.ident")

combat <- split[[1]]
comb_rna <- split[[2]]

#ChatGPT assisted with applying function across all proteomics samples
pairwise_distance_hist <- function(proteomics, rna, proteomics.subset.group = NA, 
                                   proteomics.subset.value, group.by.prot = NA, 
                                   rna.subset.group = NA, rna.subset.value, 
                                   group.by.rna = NA, stack.multi.plot = NA,
                                   distance.method = "euclidean", xlim1, xlim2) {
  
  calculate_distances <- function(rn, prot) {
    x <- rbind(prot, rn)
    return(distance(x, method = distance.method, mute.message = TRUE))
  }
  
  if (!is.na(rna.subset.group)) {
    rna <- subset(rna, subset = !!sym(rna.subset.group) == rna.subset.value)
  } 
  
  if (!is.na(proteomics.subset.group)) {
    proteomics <- subset(proteomics, subset = !!sym(proteomics.subset.group) == proteomics.subset.value)
  } 
  
  # prot_matrix <- t(proteomics@reductions$pca@cell.embeddings)
  # rna_matrix  <- t(rna@reductions$pca@cell.embeddings)
  prot_matrix <- (proteomics@assays$integrated$scale.data)
  rna_matrix  <- (rna@assays$integrated$scale.data)
  
  # iterate over all proteomics columns
  all_distances <- apply(prot_matrix, 2, function(prot_col) {
    apply(rna_matrix, 2, function(rna_col) {
      calculate_distances(rna_col, prot_col)
    })
  })
  
  # convert to long-form for plotting
  plot_data <- as.data.frame(as.table(all_distances))
  colnames(plot_data) <- c("rna_index", "prot_index", "distance")
  
  # add grouping variables if available
  if (!is.na(group.by.rna)) {
    plot_data$rna_group <- rna@meta.data[[group.by.rna]][plot_data$rna_index]
  }
  if (!is.na(group.by.prot)) {
    plot_data$prot_group <- proteomics@meta.data[[group.by.prot]][plot_data$prot_index]
  }
  if (!is.na(stack.multi.plot)) {
    plot_data$stack_group <- rna@meta.data[[stack.multi.plot]][plot_data$rna_index]
  }
  
  #MAIN PLOT
  if (!is.na(group.by.prot)) {
    plotty <- ggplot(plot_data, aes(distance, fill = prot_group)) + 
      geom_histogram(binwidth = 1) +
      labs(fill = group.by.prot)
    
    if (!is.na(stack.multi.plot)) {
      
      # Compute means per group
      mean_df <- plot_data %>%
        group_by(stack_group, prot_group) %>%
        summarise(mean_distance = mean(distance, na.rm = TRUE), .groups = "drop")
      
      splitty <- ggplot(plot_data, aes(distance, fill = stack_group)) +
        geom_histogram(binwidth = 1, alpha = 0.5, position = "identity") +
        facet_wrap(~prot_group) +
        
        # Add mean lines
        geom_vline(data = mean_df, aes(xintercept = mean_distance, color = stack_group),
                   linewidth = 1, linetype = "dashed", show.legend = FALSE) +
        
        # Add mean value labels
    
        
        theme_minimal(base_size = 14) +
        labs(
          x = "Distance",
          y = "Count",
          fill = stack.multi.plot,
          title = "Distribution of distances by group"
        ) +
        theme(
          legend.position = "top",
          strip.text = element_text(face = "bold")
        )
    
    } else {
      splitty <- ggplot(plot_data, aes(distance, fill = prot_group)) +
        geom_histogram(binwidth = 1, alpha = 0.7) +
        facet_wrap(~prot_group) +
        scale_fill_viridis_d() +
        theme_minimal() +
        labs(fill = group.by.prot)
    }
    
  } else if (!is.na(group.by.rna)) {
    plotty <- ggplot(plot_data, aes(distance, fill = rna_group)) + 
      geom_histogram(binwidth = 1) +
      labs(fill = group.by.rna)
    
    if (!is.na(stack.multi.plot)) {
      splitty <- ggplot(plot_data, aes(distance, fill = stack_group)) +
        geom_histogram(binwidth = 1, alpha = 0.7, position = "stack") +
        facet_wrap(~rna_group) +
        scale_fill_viridis_d() +
        theme_minimal() +
        labs(fill = stack.multi.plot)
    } else {
      splitty <- ggplot(plot_data, aes(distance, fill = rna_group)) +
        geom_histogram(binwidth = 1, alpha = 0.7) +
        facet_wrap(~rna_group) +
        scale_fill_viridis_d() +
        theme_minimal() +
        labs(fill = group.by.rna)
    }
  }
  plotty <- plotty + 
    xlim(xlim1, xlim2) + 
    ggtitle(glue("{distance.method} distance between proteomics ({proteomics.subset.group} {proteomics.subset.value}) and metacells"))
  
  if (!is.null(splitty)) {
    splitty <- splitty + 
      xlim(xlim1, xlim2) +
      ggtitle(glue("{distance.method} distance split by group"))
  }
  
  print(plotty)
  if (!is.null(splitty)) print(splitty)
  
  return(list(data = plot_data, plot = plotty, split_plot = splitty))
}

# #THIS ONE
# pdf(file="pca_pairwise_plots/sample_prot_euclidean.pdf", width = 16, height = 9)
# test1 <- pairwise_distance_hist(proteomics = proteomics, rna = seur_rna, proteomics.subset.group = NA,
#                                 proteomics.subset.value = NA, group.by.prot = "sample",
#                                 rna.subset.group = NA , rna.subset.value = NA,
#                                 group.by.rna = NA, distance = "euclidean",
#                                 stack.multi.plot = "sample", xlim = 0, xlim2 = 150)
# dev.off()

#THIS ONE
pdf(file="pca_pairwise_plots/combat_sample_prot_euclidean.pdf", width = 12, height = 9)
test1 <- pairwise_distance_hist(proteomics = combat, rna = comb_rna, proteomics.subset.group = NA,
                                proteomics.subset.value = NA, group.by.prot = "plot_umap_1",
                                rna.subset.group = NA , rna.subset.value = NA,
                                group.by.rna = NA, distance = "euclidean",
                                stack.multi.plot = "time", xlim1 = 10, xlim2 = 60)
dev.off()

# 
# pdf(file="pca_pairwise_plots/sample_prot_euclidean.pdf", width = 16, height = 9)
# test2 <- pairwise_distance_hist(proteomics = proteomics, rna = rna, proteomics.subset.group = NA,
#                                 proteomics.subset.value = NA, group.by.prot = "sample",
#                                 rna.subset.group = NA , rna.subset.value = NA,
#                                 group.by.rna = NA, distance = "euclidean",
#                                 stack.multi.plot = "time")
# dev.off()

# boot1 <- Mclust(as.numeric(unlist(test1$data$distance)))
# boot2 <- MclustBootstrap(boot1)
# dist_data <- test1$data
# dist_data2 <- test2$data

threshold_analysis_og <- function(dist_data, prot, threshold) {
  
  names <- names(table(dist_data$stack_group))
  print(names)
  ratios <- as.numeric((table(dist_data$stack_group))) / length(dist_data$distance)
  
  above_threshold <- dist_data %>% filter(prot_group == prot) %>% filter(distance > threshold)
  above_threshold_ratios <- as.numeric((table(above_threshold$stack_group))) / length(above_threshold$distance)
  above_threshold_ratios <- ((above_threshold_ratios / ratios) - 1) * 100
  
  above_table <- data.frame(names = names, ratios = above_threshold_ratios)
  
  below_threshold <- dist_data %>% filter(prot_group == prot) %>% filter(distance < threshold)
  below_threshold_ratios <- as.numeric((table(below_threshold$stack_group))) / length(below_threshold$distance)
  below_threshold_ratios <- ((below_threshold_ratios / ratios) - 1) * 100
  
  below_table <- data.frame(names = names, ratios = below_threshold_ratios)
  
  
  pdf(file = glue("pairwise_threshold_plots/{prot}_{threshold}.pdf"),
      width = 16, height = 9)
  
  # https://stackoverflow.com/questions/1330989/rotating-and-spacing-axis-labels-in-ggplot2
  below <- ggplot(below_table, aes(x = names, y = ratios)) + geom_bar(stat = 'identity') +
    labs(x = "sample", y = "% difference in ratio from total") + 
    ggtitle(glue("{prot} by ligand: below {threshold} | {length(below_threshold$distance)} total")) + theme(axis.text.x = element_text(angle = -90, hjust = 0))
  above <- ggplot(above_table, aes(x = names, y = ratios)) + geom_bar(stat = 'identity') +
    labs(x = "sample", y = "% difference in ratio from total") + 
    ggtitle(glue("above {threshold} | {length(above_threshold$distance)} total")) + theme(axis.text.x = element_text(angle = -90, hjust = 0))
  
  grid.arrange(below, above, nrow = 1)
  
  dev.off()
}

threshold_analysis <- function(dist_data, prot, threshold) {
  
  names <- names(table(dist_data$stack_group))
  print(names)
  ratios <- as.numeric((table(dist_data$stack_group))) / length(dist_data$distance)
  
  above_threshold <- dist_data %>% filter(prot_group == prot) %>% filter(distance > threshold)
  above_threshold_ratios <- as.numeric((table(above_threshold$stack_group))) / length(above_threshold$distance)
  above_threshold_ratios <- ((above_threshold_ratios / ratios) - 1) * 100
  
  above_table <- data.frame(names = names, ratios = above_threshold_ratios)
  
  below_threshold <- dist_data %>% filter(prot_group == prot) %>% filter(distance < threshold)
  below_threshold_ratios <- as.numeric((table(below_threshold$stack_group))) / length(below_threshold$distance)
  below_threshold_ratios <- ((below_threshold_ratios / ratios) - 1) * 100
  
  below_table <- data.frame(names = names, ratios = below_threshold_ratios)
  

  pdf(file = glue("pairwise_threshold_plots/{prot}_{threshold}.pdf"),
      width = 16, height = 9)
  
  # https://stackoverflow.com/questions/1330989/rotating-and-spacing-axis-labels-in-ggplot2
  below <- ggplot(below_table, aes(x = names, y = ratios)) + geom_bar(stat = 'identity') +
    labs(x = "sample", y = "% difference in ratio from total") + 
    ggtitle(glue("{prot} by ligand: below {threshold} | {length(below_threshold$distance)} total")) + theme(axis.text.x = element_text(angle = -90, hjust = 0))
  above <- ggplot(above_table, aes(x = names, y = ratios)) + geom_bar(stat = 'identity') +
    labs(x = "sample", y = "% difference in ratio from total") + 
    ggtitle(glue("above {threshold} | {length(above_threshold$distance)} total")) + theme(axis.text.x = element_text(angle = -90, hjust = 0))
  
  grid.arrange(below, above, nrow = 1)
  
  dev.off()
}


# threshold_analysis(dist_data, "Micro.G_Day.1_No.ligand", 117)
# 
# threshold_analysis(dist_data, "Micro.G_Day.3_with.ligand", 97.5)
# 
# threshold_analysis(dist_data, "Micro.G_Day.7_No.ligand", 84)
# 
# threshold_analysis(dist_data, "Micro.G_Day.7_with.ligand", 93)
# 
# threshold_analysis(dist_data, "Normal.G_Day.3_with.ligand", 82.5)
# 
# threshold_analysis(dist_data, "Normal.G_Day.3_No.ligand", 77)
# 
# threshold_analysis(dist_data, "Normal.G_Day.7_No.ligand", 82.5)

# manual_prot_groups <- c("Micro.G_Day.1_No.ligand", "Micro.G_Day.3_with.ligand",
#                         "Micro.G_Day.7_No.ligand", "Micro.G_Day.7_with.ligand",
#                         "Normal.G_Day.3_with.ligand", "Normal.G_Day.3_No.ligand",
#                         "Normal.G_Day.7_No.ligand")
# 
# manual_thresholds <- c(117, 97.5, 84, 93, 82.5, 77, 82.5)
# 
# two_way_table <- function(dist_data, prot_groups, thresholds) {
#   df <- rbind(prot_groups, thresholds)
#   out_matrix <- matrix(1:14, nrow = 2, dimnames = list(c("5hr", "72hr"), str_c(prot_groups, "  ", thresholds)))
#   
#   for (i in 1:ncol(df)) {
#     group_data <- dist_data %>% filter(prot_group==df[1,i])
#     below_threshold <- group_data %>% filter(distance < df[2,i])
#     
#     group_ratios <- as.numeric(table(group_data$stack_group)) / length(group_data$stack_group)
#     below_ratios <- as.numeric(table(below_threshold$stack_group)) / length(below_threshold$stack_group)
#     
#     final_ratios <- ((below_ratios/group_ratios) - 1) * 100
#     
#     out_matrix[1, i] <- final_ratios[1]
#     out_matrix[2, i] <- final_ratios[2]
#   }
#   print(out_matrix)
#   
#   test <- expand.grid(days = as.factor(c("5hr", "72hr")), groups = as.factor(str_c(prot_groups, "_", thresholds)))
#   test$Z <- expand.grid(t(out_matrix))
#   test$Z <- test$Z$Var1
#   plot2 <- ggplot(test, aes(x=days, y=groups, fill = Z)) +
#   geom_tile() + geom_text(aes(label=round(Z, 3)), color="white") +
#   scale_fill_viridis_c(, direction = 1) + coord_fixed() + 
#   guides(fill = guide_colourbar(barwidth = 0.5,
#                                 barheight = 10, title="two way table"))
#   
#   ggsave(glue("threshold_table_heatmap.png"), width = 5, height = 6)
#   return(out_matrix)
# }
# 
# tabley <- two_way_table(dist_data2, manual_prot_groups, manual_thresholds)
