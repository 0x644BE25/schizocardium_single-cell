######################################################
# PLOTTING HELPERS
#
# GOAL: This is just to make things easier to 
# replicate plotting styles from the manuscript and/or
# see how they were done.
######################################################


# ================= IMPORTS ==========================

library('Seurat')

# ================= PARAMS ===========================

# ================= METHODS ==========================

# ================= INIT DATA ========================

seu <- readRDS('./data/seurat_objects/rPCA_integrated.rds')
annot <- read.csv('./data/cluster_tissue_annotation_colors_label_position.csv',row.names=1)

# ================= COLORS ===========================

seu@misc$colors <- list()
seu@misc$colors$timepoint <- setNames(time_cols,levels(seu$timepoint))
seu@misc$colors$lifestage_bias <- setNames(lifestage_cols,c('early','other','late'))
seu@misc$colors$tissue <- sapply(levels(seu$tissue),function(x){unique(annot[annot$tissue==x,'base_color'])})
seu@misc$colors$tissue_text <- sapply(levels(seu$tissue),function(x){unique(annot[annot$tissue==x,'text_color'])})
seu@misc$colors$tissue_cluster <- setNames(annot$cluster_color,annot$cluster)
seu@misc$colors$cluster_text <- setNames(annot$text_color,annot$cluster)
seu@misc$plotting$colorsDF <- annot

# ================= PLOTTING FUNCTIONS ===============

# plot UMAPs colored by lifestage
stageUMAP <- function(seu,cells,top_cells=NULL,shuffle=TRUE) {
  library(ggplot2)
  df <- data.frame(Embeddings(seu,reduction='umap'))
  df$timepoint <- seu$timepoint
  df <- df[df$timepoint!='juvenileproboscis',]
  cellsout <- setdiff(Cells(seu),cells)
  df[cellsout,'timepoint'] <- NA
  
  pt.size <- DimPlot(seu,reduction='umap',combine=FALSE,raster=FALSE)[[1]]$layers[[1]]$aes_params$size
  if (shuffle) { cells <- sample(cells); top_cells <- sample(top_cells) }
  cell_order <- rev(unique(c(top_cells,cells,cellsout))) 
  df <- df[cell_order,]
  
  p <- ggplot(df,aes(x=UMAP_1,y=UMAP_2,color=timepoint)) + 
    geom_point(size=pt.size) + 
    scale_color_manual(values=seu@misc$colors$timepoint,na.value='#bbbbbb') +
    coord_fixed() + theme_void() + theme(legend.position='none')
  return(p)
}

# plot UMAPs colored by tissue/cluster
tissueUMAP <- function(seu,tissue=NA,labelsize=3) {
  library(ggplot2)
  if (is.na(tissue[1])) { tissue <- levels(seu$tissue) }
  df <- seu@misc$plotting$colorsDF
  seu$temp <- seu$global_clusters
  seu@meta.data[!seu$tissue %in% tissue,'temp'] <- NA
  allcells <- Cells(seu)[seu$timepoint!='juvenileproboscis']
  pt.size <- DimPlot(seu,reduction='umap',combine=FALSE,raster=FALSE)[[1]]$layers[[1]]$aes_params$size
  temp <- data.frame(cluster=seu$temp,all='ALL',x=Embeddings(seu,'umap')[,1],y=Embeddings(seu,'umap')[,2])[sample(length(Cells(seu))),]
  ggplot(temp[allcells,],aes(x=x,y=y)) + 
    geom_point(color=seu@misc$colors$na,size=pt.size) + 
    geom_point(data=temp[!is.na(temp$cluster),],aes(color=cluster),size=pt.size) +
    scale_color_manual(values=seu@misc$colors$tissue_cluster,na.value='transparent') +
    geom_label(data=df[df$tissue %in% tissue,],aes(x=x,y=y,label=cluster,fill=I(cluster_color),color=I(text_color)),size=labelsize,label.padding=unit(.1,'lines')) +
    coord_fixed() + theme_void() + theme(legend.position='none') +
    ggtitle(paste(tissue,collapse='/'))
}

# get appropriate contrasting text color
getTextColor <- function(hex) {
  ints <- strtoi(paste0('0x',c(substr(hex,2,3),substr(hex,4,5),substr(hex,6,7))))
  return({ if (mean(ints)>100) { '#000000' } else {'#ffffff' }})
}

# ================= SAVE DATA ========================
