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
                   reduction = "pca",
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

pcs <- 20
n.neighbors <- 12

supercells <- readRDS("../SuperCell/supercells/5_10_supercell.rds")

# Microarray
microarray <- read.delim('osd_324/GLDS-324_microarray_2015-10-07_Pajevic_RMA-GENE-FULL_Group.txt', )
metadata <- read.delim2('osd_324/s_GLDS-324.txt')

colnames(metadata)[2] <- "sample"
colnames(metadata)[16] <- "gravity"
colnames(metadata)[19] <- "time"

metadata$sample <- sub(" ", "_", str_c(metadata$gravity, metadata$time, sep = "_"))

microarray$Gene.Symbol <- trimws(microarray$Gene.Symbol)
colnames(microarray)[2:13] <- metadata$sample

microarray <- microarray[!microarray$Gene.Symbol == "---",]
microarray <- microarray[!is.na(microarray$Gene.Symbol),]
microarray <- microarray[!duplicated(microarray$Gene.Symbol),]

rownames(microarray) <- microarray$Gene.Symbol

micro_meta <- metadata
micro_meta <- micro_meta[,c(2,16,19)]
micro_meta$ligand <- rep("NO", 12)
micro_meta$type <- rep("Microarray", 12)

microarray <- microarray[2:13]
microarray_mat <- as.matrix(microarray)

microarray <- CreateSeuratObject(counts=microarray, assay = "microarray", meta.data = micro_meta)
microarray <- NormalizeData(microarray)
microarray <- FindVariableFeatures(microarray)
microarray <- ScaleData(microarray)
microarray <- RunPCA(microarray, npcs = 11)
microarray <- RunUMAP(microarray, dims = 1:11, n.neighbors = 5, reduction = "pca", return.model = TRUE)
saveRDS(microarray, "microarray_intersect.rds")


# Proteomics
proteomics <- readRDS("day1_3_7.rds")
proteomics <- proteomics[!grepl("_No.", colnames(proteomics))]

time <- integer(length(colnames(proteomics)))
ligand <- character(length(colnames(proteomics)))
gravity <- character(length(colnames(proteomics)))
sample <- colnames(proteomics)

time[grep("Day.1", colnames(proteomics))] <- 1
time[grep("Day.3", colnames(proteomics))] <- 3
time[grep("Day.7", colnames(proteomics))] <- 7

ligand[grep("No", colnames(proteomics))] <- "NO"
ligand[grep("with", colnames(proteomics))] <- "WITH"

gravity[grep("Micro", colnames(proteomics))] <- "MICRO"
gravity[grep("Normal", colnames(proteomics))] <- "NORMAL"

sample <- as.factor(gsub(".{2}$", "", sample))

proteomics_meta <- data.frame(time, ligand, gravity, sample)
proteomics_meta$type <- rep("Proteomics", length(proteomics))

rownames(proteomics_meta) <- colnames(proteomics)
proteomics_multi <- proteomics

# Filter out No ligand Samples
proteomics_multi <- proteomics_multi

proteomics <- CreateSeuratObject(counts=proteomics, assay = "proteomics", meta.data = proteomics_meta)

# Integration using comBat Seq

proteomics_multi <- as.data.frame(proteomics_multi)
microarray_mat <- as.data.frame(microarray_mat)

common_genes <- intersect(rownames(proteomics_multi), rownames(microarray_mat))

# Subset both data frames to only those genes
proteomics_common <- proteomics_multi[common_genes, , drop = FALSE]
microarray_common <- microarray_mat[common_genes, , drop = FALSE]

# Combine them side-by-side
merged_data <- cbind(proteomics_common, microarray_common)
merged_meta <- rbind(proteomics_meta, micro_meta)
#
merged_meta$time <- as.factor(merged_meta$time)
merged_meta$gravity <- as.factor(merged_meta$gravity)
merged_meta$sample <- as.factor(merged_meta$sample)
merged_meta$type <- as.factor(merged_meta$type)
merged_meta$ligand <- as.factor(merged_meta$ligand)

combat2 <- ComBat(merged_data, merged_meta$type)

combat <- CreateSeuratObject(counts=combat2, assay = "combat", meta.data = merged_meta)

combat@meta.data$orig.ident <- rep("ComBat", 35)
combat@meta.data$time <- (merged_meta$time)
combat@meta.data$gravity <- (merged_meta$gravity)
combat@meta.data$sample <- (merged_meta$sample)
combat@meta.data$type <- (merged_meta$type)
combat@meta.data$ligand <- (merged_meta$ligand)

intersection <- intersect(rownames(supercells$seuratMC), rownames(combat))
print(length(intersection))

combat <- subset(combat, features = intersection)
combat <- NormalizeData(combat)
combat <- FindVariableFeatures(combat)
combat <- ScaleData(combat)
combat <- RunPCA(combat, npcs = pcs)
combat <- RunUMAP(combat, dims = 1:pcs, n.neighbors = n.neighbors, reduction = "pca", return.model = TRUE)

combat$plot_1 <- combat@meta.data %>% mutate(
  composite = case_when(
    type == "combat" ~ "Cluster 1",
    ligand == "WITH" ~ "Cluster 1"
  )) %>% pull(composite)

test8 <- pca_func(combat, "sample", match_border = TRUE)
