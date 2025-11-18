# Osteocyte Multiomics
This is the code repository for my honours project, titled *Investigating the Multiomics Regulation of Space Biology for Osteocytes.* 
# Proteomics
Investigating an osteocyte proteomics dataset, using [DEP](https://bioconductor.org/packages/release/bioc/vignettes/DEP/inst/doc/DEP.html) & [Limma](https://rdrr.io/bioc/limma/) to define an Osteocyte Microgravity Gene Signature. [ClusterProfiler](https://bioconductor.org/packages//release/bioc/html/clusterProfiler.html) was used to perform functional analysis.


# Proteomics_Integration
Integrating proteomics & [publically available microarray dataset](https://osdr.nasa.gov/bio/repo/data/studies/OSD-324) together using [ComBat](https://pmc.ncbi.nlm.nih.gov/articles/PMC7518324/), a batch correction package, and scissor, 

# CellBender
Conducting quality control on an osteocyte snRNASeq dataset using [CellBender](https://github.com/broadinstitute/CellBender), a machine learning package which aims to remove the effect of ambient RNA.

# SuperCell
Analysing an osteocyte snRNASeq dataset, using [Seurat](https://satijalab.org/seurat/) as a single-cell toolbox, and using [SuperCell](https://github.com/GfellerLab/SuperCell) to reduce the sparsity of the dataset by constructing metacells.

# SuperCell Integration
Integrating metacells, proteomics, and microarray datasets together using Canonical Correlation Analysis.
