library(grid)
library(igraph)
library(Seurat)
library(glue)
library(ggplot2)
#library(ggpubr)
#library(patchwork)
#library(VennDiagram)
library(dplyr)
library(SuperCell)
library(colorspace)

#supercell <- readRDS("test/5_20_supercell.rds")

knn=5
gamma=20

# Plot Network Graphs for Time, ligand, Sample
# groups <- list(
#   time = supercell$SC$time,
#   ligand = supercell$SC$ligand,
#   sample = supercell$SC$sample)
# 
# for (group_name in names(groups)) {
#   group <- groups[[group_name]]
# 
#   q <- qualitative_hcl(length(unique(group)), palette = "Dark 3")
# 
#   png(file=glue("plots/{knn}_{gamma}_supercell_{group_name}_plot.png"), width=1080, height=1080)
#   plot5 <- supercell_plot(supercell$SC$graph.supercells,
#                           group = group,
#                           lay.method = "fr",
#                           seed = 1,
#                           main = glue("Metacells colored by {group_name}"),
#                           color.use = q
#   )
#   legend("bottomright",
#          legend = levels(as.factor(group)),  # or unique() if it's a character vector
#          col = q,
#          pch = 19,       # solid circle
#          pt.cex = 1.5,   # size of the point
#          bty = "n",      # no box around legend
#          )
#   dev.off()
# }

# Plot network graphs for sample, ligand, stiffness
for (day in list(c(hr_5, hr_72))) {
  hours <- unique(day$seuratMC$time)
  groups <- list(
    stiffness = day$SC$stiffness,
    ligand = day$SC$ligand,
    sample = day$SC$sample)
  for (group_name in names(groups)) {
    group <- groups[[group_name]]
    
    q <- qualitative_hcl(length(unique(group)), palette = "Dark 3")
    
    png(file=glue("plot_day/{hours}_supercell_{group_name}_plot.png"), width=1080, height=1080)
    plot5 <- supercell_plot(day$SC$graph.supercells,
                            group = group,
                            lay.method = "fr ",
                            seed = 1,
                            main = glue("{hours} hours Metacells colored by {group_name}"),
                            color.use = q
    )
    legend("bottomright",
           legend = levels(as.factor(group)),  # or unique() if it's a character vector
           col = q,
           pch = 19,       # solid circle
           pt.cex = 1.5,   # size of the point
           bty = "n",      # no box around legend
    )
    dev.off()
  }
}

# Plot network graphs w/ different methods
# for (method in c("nicely", "fr", "components", "drl", "graphopt")) {
#   
#   group = supercell$SC$sample
#   q <- qualitative_hcl(length(unique(group)), palette = "Dark 3")
#   
#   png(file=glue("plot_network_methods/{knn}_{gamma}_supercell_{method}_plot.png"), width=1080, height=1080)
#   plot5 <- supercell_plot(supercell$SC$graph.supercells,
#                           group = group,
#                           lay.method = method,
#                           seed = 1,
#                           
#                           main = glue("{method}"),
#                           color.use = q
#   )
#   legend("bottomright",
#          legend = levels(as.factor(group)),  # or unique() if it's a character vector
#          col = q,
#          pch = 19,       # solid circle
#          pt.cex = 1.5,   # size of the point
#          bty = "n",      # no box around legend
#   )
#   dev.off()
# }




