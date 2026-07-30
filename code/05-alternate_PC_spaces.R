######################################################
# ALTERNATE PC SPACES
#
# GOAL: USe protein domains to identify TF and non-TF
# gene sets.
# 
# Use rbh BLAST against mouse matches + gprofiler 
# mouse annotations to find "metabolic" and 
# "catabolic" genes sets.
#
# Once gene sets are determined,do PCA and generate 
# UMAP embeddings each gene set to asses cell type and
# life stage separation.
#
# GO terms list downloaded 2026.01.30@14:45PST from
# https://biit.cs.ut.ee/gprofiler/gost
# detailed version info available in 
# ./data/go_terms/gprofiler_full_mmusculus_version_info.txt
######################################################

if (!dir.exists('./data/alternate_PCs')) { dir.create('./data/alternate_PCs') }

# ================= IMPORTS ==========================

library(Seurat)
library(SeuratWrappers)
library(gprofiler2)
library(ggplot2)

# ================= PARAMS ===========================

maxPCs <- 55

stageCols <- c('earlylarva'='#FFCC22','latelarva'='#C09000',
               'shared'='#CCCCCC','duringmeta'='#666666',
               'earlyjuvenile'='#66BBFF','latejuvenile'='#3075c0')

# ================= METHODS ==========================

# ================= INIT DATA ========================

seu <- readRDS('./data/seurat_objects/rPCA_integrated.rds')
seu <- subset(seu,cells=Cells(seu)[seu$timepoint!='juvenileproboscis'])
key <- read.csv('./data/SCA_to_geneName.csv',row.names=1)

# TFs and non-TFs
tfs <- read.delim('./data/transcription_factor_domains.tsv',sep='\t',header=FALSE)[,1]


# meta- cata- bolic genes
blast <- data.table::fread('./data/go_terms/mmus_reciprocal_best_hits.csv.gz')
curated_terms <- read.delim('./data/go_terms/GO_terms_from_Paul_20240919.txt',header=TRUE,sep=';',row.names=1)
# NOTE: you'll need to download the following yourself, see
# ./data/go_terms/gprofiler_full_mmusculus_version_info.txt
# for a full list of params
terms <- readLines('./data/go_terms/gprofiler_full_mmusculus.ENSG.gmt')

# ================= GET GENES ========================

# TFs and non-TFs
tfs <- unique(sapply(strsplit(tfs,'\\.'),function(x){ x[1] }))
tfs <- intersect(rownames(key)[key$SCA %in% tfs],rownames(seu[['integrated']]))
nontfs <- setdiff(rownames(seu[['integrated']]),tfs)
writeLines(tfs,'./data/alternate_PCs/TF_genes.txt')
writeLines(nontfs,'./data/alternate_PCs/non-TF_genes.txt')

# meta- and cata- bolic genes
ensg <- gconvert(unique(substr(blast$mouse,1,18)),organism='mmusculus',target='ENSG')
rownames(ensg) <- ensg$input
blast$gene_id <- ensg[substr(blast$mouse,1,18),'target']
terms <- do.call(rbind,lapply(strsplit(terms,'\t'),function(x) {
  data.frame(term=x[1],description=x[2],genes=paste(x[3:length(x)],collapse=','))
}))
rownames(terms) <- terms$term
terms <- terms[startsWith(terms$term,'GO:'),]

ab_terms <- terms[grepl('aboli',terms$description,ignore.case=TRUE),1:2]
write.csv(ab_terms,'./data/go_terms/mmus_meta_cata_bolic_GO_terms.csv',row.names=FALSE)

meta_terms <- ab_terms[grepl('metaboli',ab_terms$description,ignore.case=TRUE),]
cata_terms <- ab_terms[grepl('cataboli',ab_terms$description,ignore.case=TRUE),]

meta_genes <- unique(unlist(lapply(meta_terms$term,function(term){
  mmus_genes <- strsplit(terms[term,'genes'],',')[[1]]
  return(blast[blast$gene_id %in% mmus_genes,'SCA'])
})))
meta_genes <- intersect(rownames(key)[key$SCA %in% meta_genes],rownames(seu[['integrated']]))

cata_genes <- setdiff(unlist(lapply(cata_terms$term,function(term){
  mmus_genes <- strsplit(terms[term,'genes'],',')[[1]]
  return(blast[blast$gene_id %in% mmus_genes,'SCA'])
})),tfs)
cata_genes <- intersect(rownames(key)[key$SCA %in% cata_genes],rownames(seu[['integrated']]))

catgo <- curated_terms[curated_terms$term_name=='Catabolic process','term_id']
catgo_genes <- unlist(strsplit(terms[catgo,'genes'],',')[[1]])
catgo_genes <- unique(blast[blast$gene_id %in% unique(unlist(strsplit(terms[catgo,'genes'],',')[[1]])),'SCA'])
catgo_genes <- intersect(rownames(key)[key$SCA %in% catgo_genes],rownames(seu[['integrated']]))

writeLines(meta_genes,'./data/alternate_PCs/metabolic_genes.txt')
writeLines(cata_genes,'./data/alternate_PCs/catabolic_genes.txt')
writeLines(catgo_genes,paste0('./data/alternate_PCs/catabolic_process_',catgo,'_genes.txt'))

# ================= PCA ==============================

tf_genes <- readLines('./data/alternate_PCs/TF_genes.txt')
nontf_genes <- readLines('./data/alternate_PCs/non-TF_genes.txt')
meta_genes <- readLines('./data/alternate_PCs/metabolic_genes.txt')
cata_genes <- readLines('./data/alternate_PCs/catabolic_genes.txt')
catgo_genes <- readLines('./data/alternate_PCs/catabolic_process_GO:0009056_genes.txt')
geneList <- list(
  'TF'=tf_genes,
  'non-TF'=nontf_genes,
  'catabolic_GO0009056'=catgo_genes,
  'metabolic_and_catabolic'=c(meta_genes,cata_genes),
  'non-metabolic_non-catabolic'=setdiff(rownames(seu[['integrated']]),c(meta_genes,cata_genes)))

int <- seu[['integrated']]

for (n in names(geneList)) {
  genes <- geneList[[n]]
  cat(n,genes[1:5],'\n')
  hm <- CreateAssayObject(counts=int$data[genes,])
  hm@meta.features <- int@meta.features[genes,]
  hm$data <- hm$counts
  curr <- seu
  curr[['temp']] <- hm
  curr@active.assay <- 'temp'
  curr <- FindVariableFeatures(curr,assay='temp',nfeatures=min(2000,length(genes)))
  curr <- RunPCA(ScaleData(curr),npcs=maxPCs)
  curr <- FindNeighbors(curr,dims=1:maxPCs,graph.name=c('curr_NN','curr_SNN'))
  curr <- RunUMAP(curr,dims=1:maxPCs,reduction='pca',reduction.name='xxx',reduction.key='xxx_')
  p <- DimPlot(curr,reduction='xxx',group.by='timepoint',raster=FALSE) +
    scale_color_manual(values=stageCols) +
    theme_void() + coord_fixed() + ggtitle(n)
  ggsave(paste0('./figures/UMAPs_',n,'_timepoint.png'),p,width=10,height=8,units='in',dpi=350)
  write.csv(Embeddings(curr,reduction='pca'),paste0('./data/alternate_PCs/',n,'_55PC_cell_embeddings.csv'),row.names=TRUE)
  write.csv(Embeddings(curr,reduction='xxx'),paste0('./data/alternate_PCs/',n,'_UMAP_cell_embeddings.csv'),row.names=TRUE)
  write.csv(curr[['pca']]@feature.loadings,paste0('./data/alternate_PCs/',n,'_55PC_feature_loadings.csv'),row.names=TRUE)
}

