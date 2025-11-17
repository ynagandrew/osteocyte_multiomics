library(glue)
library(Seurat)
library(ggplot2)
library(gridExtra)
library(dplyr)

#supercells <- readRDS("C:/Users/Andrew Yang/Desktop/proteomics/multiomics/5_10_supercell_w_scissor.rds")

names <- list(c("ligand_scissor", "No Ligand", "Ligand"),
           c("gravity_scissor", "Normal Gravity", "Microgravity"),
           c("time1_scissor", "Day-Dissimilar", "Day-Similar"),
           c("time3_scissor", "Day3-", "Day3+"),
           c("time7_scissor", "Day7-", "Day7+"))

for (name in names) {
  neglabel <- name[[2]]
  poslabel <- name[[3]]
  name <- name[[1]]

  pdf(file = glue("thesis_scissor_plots/{name}.pdf"), width=10, height=10)
  sc <- DimPlot(supercells, reduction = 'umap', group.by = name, cols = c('grey','indianred1','green'), pt.size = 1.2, order = c(2,1))

  df <- supercells@meta.data %>%
    select(sample, !!sym(name)) %>%
    mutate(!!sym(name) := as.factor(!!sym(name))) %>%  # 🔹 convert to factor
    group_by(sample, !!sym(name)) %>%
    summarise(count = n(), .groups = "drop") %>%
    group_by(sample) %>%
    mutate(prop = count / sum(count))

  print(df)

  # stacked bar plot
  sample_dist <- ggplot(df, aes(x = sample, y = prop, fill = !!sym(name))) +
    geom_bar(stat = "identity", position = "stack", color = "black") +
    scale_fill_manual(values = c('grey', 'indianred1', 'green'),
                      labels = c("No Phenotype", neglabel, poslabel)) +
    labs(
      x = "sample",
      y = "Proportion of cells",
      fill = name
    ) +
    theme_classic(base_size = 14) +
    theme(
      legend.position = "right",
      axis.text.x = element_text(angle = 45, hjust = 1),
      axis.title.x.top = element_text(size = 14, face = "bold")
    )

  df <- supercells@meta.data %>%
    select(ligand, !!sym(name)) %>%
    mutate(!!sym(name) := as.factor(!!sym(name))) %>%  # 🔹 convert to factor
    group_by(ligand, !!sym(name)) %>%
    summarise(count = n(), .groups = "drop") %>%
    group_by(ligand) %>%
    mutate(prop = count / sum(count))

  print(df)

  # stacked bar plot
  ligand_dist <- ggplot(df, aes(x = ligand, y = prop, fill = !!sym(name))) +
    geom_bar(stat = "identity", position = "stack", color = "black") +
    scale_fill_manual(values = c('grey', 'indianred1', 'green')) +
    labs(
      x = "ligand",
      y = "Proportion of cells",
      fill = name
    ) +
    theme_classic(base_size = 14) +
    theme(
      legend.position = "none",
      axis.text.x = element_text(angle = 45, hjust = 1),
      axis.title.x.top = element_text(size = 14, face = "bold")
    )

  df <- supercells@meta.data %>%
    select(stiffness, !!sym(name)) %>%
    mutate(!!sym(name) := as.factor(!!sym(name))) %>%  # 🔹 convert to factor
    group_by(stiffness, !!sym(name)) %>%
    summarise(count = n(), .groups = "drop") %>%
    group_by(stiffness) %>%
    mutate(prop = count / sum(count))

  print(df)

  # stacked bar plot
  stiffness_dist <- ggplot(df, aes(x = stiffness, y = prop, fill = !!sym(name))) +
    geom_bar(stat = "identity", position = "stack", color = "black") +
    scale_fill_manual(values = c('grey', 'indianred1', 'green')) +
    labs(
      x = "stiffness",
      y = "Proportion of cells",
      fill = name
    ) +
    theme_classic(base_size = 14) +
    theme(
      legend.position = "none",
      axis.text.x = element_text(angle = 45, hjust = 1),
      axis.title.x.top = element_text(size = 14, face = "bold")
    )

  df <- supercells@meta.data %>%
    select(time, !!sym(name)) %>%
    mutate(!!sym(name) := as.factor(!!sym(name))) %>%  # 🔹 convert to factor
    group_by(time, !!sym(name)) %>%
    summarise(count = n(), .groups = "drop") %>%
    group_by(time) %>%
    mutate(prop = count / sum(count))

  print(df)

  # stacked bar plot
  time_dist <- ggplot(df, aes(x = time, y = prop, fill = !!sym(name))) +
    geom_bar(stat = "identity", position = "stack", color = "black") +
    scale_fill_manual(values = c('grey', 'indianred1', 'green')) +
    labs(
      x = "time",
      y = "Proportion of cells"
    ) +
    theme_classic(base_size = 14) +
    theme(
      legend.position = "none",
      axis.text.x = element_text(angle = 45, hjust = 1),
      axis.title.x.top = element_text(size = 14, face = "bold")
    )

  groblist <- list(sc, sample_dist, ligand_dist, stiffness_dist, time_dist)

  grid.arrange(
    grobs = groblist,
    layout_matrix = rbind(c(1,1,1,1,3),
                          c(1,1,1,1,4),
                          c(1,1,1,1,5),
                          c(2,2,2,2,2))
  )

  dev.off()


}

for (name in names) {
  neglabel <- name[[2]]
  poslabel <- name[[3]]
  name <- name[[1]]
  
  df <- supercells@meta.data %>%
        select(time, !!sym(name)) %>%
        mutate(!!sym(name) := as.factor(!!sym(name))) %>%  # 🔹 convert to factor
        group_by(time, !!sym(name)) %>%
        summarise(count = n(), .groups = "drop") %>%
        group_by(time) %>%
        mutate(prop = count / sum(count))

      print(df)
  
  # stacked bar plot
      time_dist <- ggplot(df, aes(x = time, y = prop, fill = !!sym(name))) +
        geom_bar(stat = "identity", position = "stack", color = "black") +
        
        # Add proportion labels
        geom_text(
          aes(label = ifelse(prop > 0.03, scales::percent(prop, accuracy = 1), "")),
          position = position_stack(vjust = 0.5),
          size = 3.5
        ) +
        
        scale_fill_manual(values = c('grey', 'indianred1', 'green')) +
        labs(
          x = "time",
          y = "Proportion of cells"
        ) +
        theme_classic(base_size = 14) +
        theme(
          legend.position = "none",
          axis.text.x = element_text(angle = 45, hjust = 1),
          axis.title.x.top = element_text(size = 14, face = "bold")
        )
      
  
  groblist <- list(sc, sample_dist, ligand_dist, stiffness_dist, time_dist)
  
  pdf(file=glue("thesis_scissor_plots/{name}_bar.pdf"))
  grid.arrange(time_dist)
  
  dev.off()
  
  
}

for (name in names) {
  supercells <- SetIdent(supercells, value=supercells@meta.data[[name]])
  test <- FindMarkers(object = supercells, ident.1 = 2, ident.2 = 1)
  print(head(test))
}



