######################################################
# GO TERM ANALYSIS
#
# GOAL: Find differentially-expressed genes at various
# levels of granularity, convert to mmus genes using
# BLAST reciprocal best hit matches, then use for 
# GO term enrichment analysis with gprofiler.
######################################################

setwd('/n/projects/cb2350/schizocardium_sc/publication_pipeline/')

# ================= IMPORTS ==========================

library(Seurat)
library(gprofiler2)

# ================= PARAMS ===========================

maxP <- 10^(-4)

# ================= METHODS ==========================

getMarkers <- function(seu,meta,genes=NULL) {
  seu@active.ident <- meta
  cat('\n')
  ml <- lapply(levels(meta),function(x){
    cat(' ',x)
    m <- FindMarkers(seu,ident.1=x,only.pos=TRUE,assay='SCT',features=genes)
    return(m[m$p_val_adj<maxP & m$avg_log2FC>=log2(1.25),])
  })
  names(ml) <- levels(meta)
  cat('\n')
  return(ml)
}


# ================= INIT DATA ========================

seu <- readRDS('../data/seurat_objects/rPCA_integrated.rds')
seu <- subset(seu,cells=Cells(seu)[seu$timepoint!='juvenileproboscis'])
seu$timepoint <- droplevels(seu$timepoint)

blast <- read.csv('./data/go_terms/mmus_reciprocal_best_hits.csv')
annot <- read.delim('./data/SCA_alias_name.txt',header=FALSE,col.names=c('sca','gene','desc'))
rownames(annot) <- annot$sca

terms <- read.delim('./data/go_terms/GO_terms_from_Paul_20240919.txt',sep=';')

# ================= GENERATE GO DATA =================

# TIMEPOINT ======================
tpMarkers <- getMarkers(seu,meta=seu$timepoint,genes=NULL)
tpGenes <- setNames(lapply(tpMarkers,function(m){
  m <- m[order(m$avg_log2FC,decreasing=TRUE),]
  do.call(rbind,lapply(rownames(m),function(x){
    mm <- substr(unlist(blast[blast$gene==x,'mouse']),1,18)
    if (length(mm)>0) { data.frame(schizo=x,mouse=mm) } else { NULL }
  }))
}),names(tpMarkers))
tpMus <- sapply(tpGenes,function(x){ x$mouse })
saveRDS(gost(tpMus,organism='mmusculus',ordered_query=TRUE,multi_query=TRUE),'./data/go_terms/GOST_results_timepoint.Rds')

# TISSUE =========================
tisMarkers <- getMarkers(seu,meta=seu$tissue,genes=NULL)
tisGenes <- setNames(lapply(tisMarkers,function(m){
  m <- m[order(m$avg_log2FC,decreasing=TRUE),]
  do.call(rbind,lapply(rownames(m),function(x){
    mm <- substr(unlist(blast[blast$gene==x,'mouse']),1,18)
    if (length(mm)>0) { data.frame(schizo=x,mouse=mm) } else { NULL }
  }))
}),names(tisMarkers))
tisMus <- sapply(tisGenes,function(x){ x$mouse })
saveRDS(gost(tisMus,organism='mmusculus',ordered_query=TRUE,multi_query=TRUE),'./data/go_terms/GOST_results_tissue.Rds')

# GLOBAL CLUSTER =================
gcMarkers <- getMarkers(seu,meta=seu$global_clusters,genes=NULL)
gcGenes <- setNames(lapply(gcMarkers,function(m){
  m <- m[order(m$avg_log2FC,decreasing=TRUE),]
  do.call(rbind,lapply(rownames(m),function(x){
    mm <- substr(unlist(blast[blast$gene==x,'mouse']),1,18)
    if (length(mm)>0) { data.frame(schizo=x,mouse=mm) } else { NULL }
  }))
}),names(gcMarkers))
gcMus <- sapply(gcGenes,function(x){ x$mouse })
saveRDS(gost(gcMus,organism='mmusculus',ordered_query=TRUE,multi_query=TRUE),'./data/go_terms/GOST_results_global_cluster.Rds')
