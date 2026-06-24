######################################################
# CLUSTER MARKERS AND ANNOTATIONS
#
# GOALS: Find markers for all clusters.
#
# Add cluster cell type annotations, along with colors
# for plotting.
#
# Determine life stage bias for each cluster and add 
# to metadata, along with colors for plotting.
######################################################

# ================= IMPORTS ==========================

library(Seurat)
library(ggplot2)
library(openxlsx)

# ================= PARAMS ===========================

stages <- c('EL'='earlylarva','LL'='latelarva',
             'DM'='duringmeta',
             'EJ'='earlyjuvenile','LJ'='latejuvenile',
             'JP'='juvenileproboscis')
stages_include <- setdiff(stages,'juvenileproboscis')

# ================= METHODS ==========================


# ================= INIT DATA ========================

seu <- readRDS('./data/seurat_objects/rPCA_integrated.rds')
annot <- read.csv('./data/cluster_tissue_annotation_colors_label_position.csv',row.names=1)

# ================= CLUSTER MARKERS ==================

seu@active.ident <- seu$global_clusters
m <- FindAllMarkers(m)
write.csv(m,'/data/markers/rPCA_integrated_global_cluster_markers.csv')

# ================= CELL TYPE ========================

seu$tissue <- factor(annot[as.character(seu$global_clusters),'tissue'],levels=sort(unique(annot$tissue)))

# ================= LIFESTAGE ========================

seu$timepoint <- setNames(factor(stages[substr(Cells(seu),1,2)],levels=stages),Cells(seu))

lsb <- do.call(rbind,lapply(levels(seu$global_clusters),function(cl){
  cells <- Cells(seu)[seu$global_clusters==cl]
  curr <- table(seu$timepoint[cells])[stages_include]/table(seu$timepoint)[stages_include]
  curr <- curr/sum(curr)
  top <- c(early=sum(curr[c('earlylarva','latelarva')]),
           intermediate=as.numeric(curr['duringmeta']),
           late=sum(curr[c('earlyjuvenile','latejuvenile')]))
  top <- names(top)[which.max(top)]
  df <- data.frame(cluster=cl,top=top)
  for (x in names(curr)) { df[,x] <- curr[x] }
  return(df)
}))
lsb[lsb$cluster=='32','top'] <- 'intermediate'
write.csv(lsb,'./data/cluster_livestage_composition.csv',row.names=FALSE)

seu$lifestage_bias <- lsb[as.character(seu$global_clusters),'top']
seu$lifestage_bias <- factor(seu$lifestage_bias,levels=c('early','intermediate','late'))

# ================= SAVE RDS =========================

saveRDS(seu,'./data/seurat_objects/rPCA_integrated.rds',compress=FALSE)
