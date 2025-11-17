library(devtools)
library(CWGCNA)
library(doParallel)
registerDoParallel(cores=8)

# allost <- readRDS("../SuperCell/Osteocyte_v2025.rds")
supercells <- readRDS("../SuperCell/supercells/5_10_supercell.rds")

# Microarray
microarray <- read.delim('osd_324/GLDS-324_microarray_2015-10-07_Pajevic_RMA-GENE-FULL_Group.txt', )
metadata <- read.delim2('osd_324/s_GLDS-324.txt')

colnames(metadata)[2] <- "sampleid"
colnames(metadata)[16] <- "gravity"
colnames(metadata)[19] <- "time"

microarray$Gene.Symbol <- trimws(microarray$Gene.Symbol)
colnames(microarray)[2:13] <- metadata$sample

microarray <- microarray[!microarray$Gene.Symbol == "---",]
microarray <- microarray[!is.na(microarray$Gene.Symbol),]
microarray <- microarray[!duplicated(microarray$Gene.Symbol),]

rownames(microarray) <- microarray$Gene.Symbol

microarray <- microarray[2:13]
microarray <- as.matrix(microarray)

metadata <- metadata[,c(2,16,19)]

# Proteomics
proteomics <- readRDS("day1_3_7.rds")
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

time <- as.factor(time)
ligand <- as.factor(ligand)
gravity <- as.factor(gravity)

meta.data <- data.frame(time, ligand, gravity, sample)

rownames(meta.data) <- colnames(proteomics)

microarray <- featuresampling(betas = microarray,
                            topfeatures = 10000,
                            variancetype = "sd",
                            threads = 8)

anovares <- featuresampling(betas = microarray$betas,
                            pddat = metadata,
                            anova = TRUE,

                            plotannovares = TRUE,
                            featuretype = "probe",
                            plottitlesuffix = "placenta",
                            titlesize = 18,
                            textsize = 16,
                            threads = 8)

testing <- diffwgcna(dat = microarray$betas, 
                     pddat = metadata, 
                     responsevarname = "gravity", 
                     confoundings = c("time"), 
                     featuretype = "gene", 
                     
                     topvaricancetype = "sd", 
                     topvaricance = 5000, 
                     
                     powers = seq(1, 20, 1), 
                     rsqcutline = 0.8, 
                     
                     mediation = TRUE, 
                     
                     topn = 1, 
                     
                     plot = FALSE, 
                     labelnum = 5, 
                     titlesize = 17, 
                     textsize = 16, 
                     annotextsize = 5, 
                     
                     pvalcolname = "adj.P.Val", 
                     pvalcutoff = 0.05, 
                     isbetaval = TRUE, 
                     absxcutoff = 0, 
                     diffcutoff = 0, 
                     
                     threads = 8)