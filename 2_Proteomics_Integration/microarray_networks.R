

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