library(clusterProfiler)
library(org.Mm.eg.db)
library(msigdbr)
library(enrichplot)
library(stringr)
library(gridExtra)
library(ggplot2)
library(ReactomePA)

os_micro_sig_name <- readRDS("osteocyte_microgravity_signature_name.rds")
os_micro_sig_uniprot <- readRDS("osteocyte_microgravity_signature_uniprot.rds")


universe_name <- readRDS("universe_name.rds")
universe_uniprot <- readRDS("universe_uniprot.rds")

universe_name <- str_split_i(universe_name, ";", 1)

# os_micro_sig_entrez <- bitr(os_micro_sig_name, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Mm.eg.db")$ENTREZID
# universe_entrez <- bitr(universe_name, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Mm.eg.db")$ENTREZID
# 
# mm_msigdb_df <- msigdbr(species = "Mus musculus", collection = 'H')

kegg <- enrichKEGG(gene = os_micro_sig_uniprot,
                   organism = "mmu",
                   keyType = "uniprot",
                   pvalueCutoff = 0.05,
                   universe = universe_uniprot,
)

cc <- enrichGO(gene = os_micro_sig_uniprot,
                 universe = universe_uniprot,
                 OrgDb = org.Mm.eg.db,
                 ont = "CC",
                 pvalueCutoff = 0.05,
                 pAdjustMethod = "BH",
                 keyType = "UNIPROT"
                 )

bp <- enrichGO(gene = os_micro_sig_uniprot,
               universe = universe_uniprot,
               OrgDb = org.Mm.eg.db,
               ont = "BP",
               pvalueCutoff = 0.05,
               pAdjustMethod = "BH",
               keyType = "UNIPROT"
)

mf <- enrichGO(gene = os_micro_sig_uniprot,
               universe = universe_uniprot,
               OrgDb = org.Mm.eg.db,
               ont = "MF",
               pvalueCutoff = 0.05,
               pAdjustMethod = "BH",
               keyType = "UNIPROT"
)

lig_bp <- enrichGO(gene = ligand_list,
               universe = universe_uniprot,
               OrgDb = org.Mm.eg.db,
               ont = "BP",
               pvalueCutoff = 0.05,
               pAdjustMethod = "BH",
               keyType = "UNIPROT"
)

lig_cc <- enrichGO(gene = ligand_list,
                   universe = universe_uniprot,
                   OrgDb = org.Mm.eg.db,
                   ont = "CC",
                   pvalueCutoff = 0.05,
                   pAdjustMethod = "BH",
                   keyType = "UNIPROT"
)

# hallmark <- enricher(gene = os_micro_sig_name,
#                      universe = universe_name,
#                      pvalueCutoff = 0.05,
#                      pAdjustMethod = "BH",
#                      TERM2GENE = dplyr::select(mm_msigdb_df, gs_name, gene_symbol)
# 
# )
# 
# reactome <- enrichPathway(gene = os_micro_sig_entrez, 
#                           organism = "mouse", 
#                           universe = universe_entrez, 
#                           pvalueCutoff = 0.05,
#                           pAdjustMethod = "BH",
#                           readable = TRUE)

pdf(file="ora_plots/cc.pdf", width=16, height=9)
grid.arrange(dotplot(cc), upsetplot(cc), nrow=1, top="cc ORA")
dev.off()

pdf(file="ora_plots/bp.pdf", width=16, height=9)
grid.arrange(dotplot(bp), upsetplot(bp), nrow=1, top="bp ORA")
dev.off()

pdf(file="ora_plots/mf.pdf", width=16, height=9)
grid.arrange(dotplot(mf), upsetplot(mf), nrow=1, top="mf ORA")
dev.off()

pdf(file="ora_plots/kegg.pdf", width=16, height=9)
grid.arrange(dotplot(kegg), upsetplot(kegg), nrow=1, top="kegg ORA")
dev.off()

# pdf(file="ora_plots/hallmark.pdf", width=16, height=9)
# grid.arrange(dotplot(hallmark), upsetplot(hallmark), nrow=1, top="kegg ORA")
# dev.off()
# 
# pdf(file="ora_plots/reactome.pdf", width=16, height=9)
# grid.arrange(dotplot(reactome), upsetplot(reactome), nrow=1, top="reactome ORA")
# dev.off()

# CC Notable Terms: Postsynapse, myelin sheath,
# Cytosolic ribosome, cytosolic large ribosomal subunit, cytosolic small ribosomal subunit
# Cytoplasmic stress granule
cnetplot(cc)

# BP Notable Terms: Regulation of cellular localisation, cytoplasmic translation, response to mechanical stimulus
# translation at synapse, translation at presynapse, translation at postsynapse
cnetplot(mf)

# MF Notable Terms: hydrolase activity, acting on acid anhydrides ???

# Kegg Notable Terms: None

# Hallmark Notable Terms: None

# Reactome Notable Terms: Nonsense-mediated decay?



