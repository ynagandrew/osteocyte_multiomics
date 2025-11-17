library(Seurat)
library(ggplot2)
library(grid)
library(gridExtra)
library(MOFA2)
library(data.table)
library(tidyr)
library(dplyr)
library(reticulate)
library(stringr)

# Multiomics Factor Analysis Seurat
# seurat.combined <- readRDS(file="seurat.combined.rds")
# mofa_input <- seurat.combined[[1]]
# DefaultAssay(mofa_input) <- "proteomics"
# mofa_input[["RNA"]] <- as(object = mofa_input[["RNA"]], Class = "Assay")
# mofa_input[["proteomics"]] <- as(object = mofa_input[["proteomics"]], Class = "Assay")
# mofa_input[["integrated"]] <- NULL
# 
# test1 <- create_mofa_from_Seurat(mofa_input, assays = c("RNA", "proteomics"))

# proteomics <- readRDS("proteomics_reduct.rds")
# supercells <- readRDS("../SuperCell/supercells/5_10_supercell.rds")
# supercells <- supercells$seuratMC
# 
# merged <- merge(supercells, y = proteomics, project = "merged")
# 
# merged <- NormalizeData(merged, normalization.method = "LogNormalize", assay = "proteomics")
# merged <- ScaleData(merged,assay="proteomics")
# merged <- FindVariableFeatures(merged, assay = "proteomics")
# 
# merged <- NormalizeData(merged, normalization.method = "LogNormalize", assay = "RNA")
# merged <- ScaleData(merged,assay="RNA")
# merged <- FindVariableFeatures(merged, assay = "RNA")
# 
# merged[["RNA"]] <- as(object = merged[["RNA"]], Class = "Assay")
# merged[["proteomics"]] <- as(object = merged[["proteomics"]], Class = "Assay")

# prot_vars <- merged@assays[["proteomics"]]@meta.data[["var.features"]]
# prot_vars <- prot_vars[!is.na(prot_vars)]
# 
# rna_vars <- merged@assays[["RNA"]]@meta.data[["var.features"]]
# rna_vars <- rna_vars[!is.na(rna_vars)]
# 
# variable_features = list("RNA" = rna_vars, 
#                       "proteomics" = prot_vars)
# mofa3 <- create_mofa_from_Seurat(merged, features = variable_features)

# Multiomics Factor Analysis Matrix
# supercells <- readRDS("../SuperCell/supercells/5_10_supercell.rds")
# proteomics <- readRDS("day1_3_7.rds")
# 
# prot_mat <- as.matrix(proteomics)
# supe_mat <- as.matrix(supercells$seuratMC@assays$RNA@layers$scale.data)
# 
# data = c(prot_mat, supe_mat)
# 
# print(Sys.time())
# mofa1 <- create_mofa_from_matrix(data)

# Multiomics Factor Analysis Data.frame

proteomics <- readRDS("proteomics_reduct.rds")
supercells <- readRDS("../SuperCell/supercells/5_10_supercell.rds")
supercells <- supercells$seuratMC

merged <- merge(supercells, y = proteomics, project = "merged")

merged <- NormalizeData(merged, normalization.method = "LogNormalize", assay = "proteomics")
merged <- ScaleData(merged,assay="proteomics")
merged <- FindVariableFeatures(merged, assay = "proteomics")

merged <- NormalizeData(merged, normalization.method = "LogNormalize", assay = "RNA")
merged <- ScaleData(merged,assay="RNA")
merged <- FindVariableFeatures(merged, assay = "RNA")

rna_values <- merged@assays[["RNA"]]@layers[["scale.data"]]   # genes × cells
rna_features <- rownames(rna_values)
rna_view <- "RNA"

# Transpose so samples are rows, features are columns
df <- data.frame(t(rna_values))

rdf <- df %>%
  tibble::rownames_to_column(var = "sample") %>%
  pivot_longer(
    cols = -sample,                     # pivot all features
    names_to = "feature",
    values_to = "value"
  ) %>%
  mutate(view = rna_view,
         sample = paste(sample, supercells$sample[match(sample, names(supercells$sample))])) %>%
  select(sample, feature, view, value)

head(rdf)

prot_values <- merged@assays[["proteomics"]]@layers[["scale.data"]] 
prot_features <- rownames(prot_values)
prot_view <- "proteomics"

# Transpose so samples are rows, features are columns
df <- data.frame(t(prot_values))

pdf <- df %>%
  tibble::rownames_to_column(var = "sample") %>%
  pivot_longer(
    cols = -sample,                     # pivot all features
    names_to = "feature",
    values_to = "value"
  ) %>%
  mutate(view = prot_view,
         sample = paste(sample, proteomics[[1]]$sample[match(sample, names(proteomics[[1]]$sample))])) %>%
  select(sample, feature, view, value)

head(rdf)

mofa_input <- bind_rows(pdf, rdf)
mofa1 <- create_mofa_from_df(mofa_input)

data_opts <- get_default_data_options(mofa1)
data_opts$scale_views <- TRUE
head(data_opts)

model_opts <- get_default_model_options(mofa1)
head(model_opts)

train_opts <- get_default_training_options(mofa1)
head(train_opts)

MOFAobject <- prepare_mofa(
  object = mofa1,
  data_options = data_opts, 
  model_options = model_opts,
  training_options = train_opts
)

outfile = file.path(getwd(),"model.hdf5")
MOFAobject.trained <- run_mofa(MOFAobject, outfile, use_basilisk = TRUE)

new_sample <- apply(X = MOFAobject.trained@samples_metadata, MARGIN = 1, 
                    FUN = str_split_i, pattern = " ", i = 2)
new_sample <- new_sample[2,]

split_sample <- apply(X = new_sample, MARGIN = 1, FUN = str_split, pattern = "_")
split_sample <- split_sample[2]
