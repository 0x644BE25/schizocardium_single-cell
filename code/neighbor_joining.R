######################################################
# NEIGHBOR JOINING
#
# GOAL: Use cluster centroids in PC space to perform
# neighbor joining.
######################################################

# ================= IMPORTS ==========================

library(Seurat)
library(ggplot2)
library(ape)
# BiocManager::install("ggtree")
library(ggtree)

# ================= PARAMS ===========================

# ================= METHODS ==========================

getTissue <- function(x) {
  tis <- sapply(x,function(xx){
    if (is.na(xx)) { return(xx) }
    return(as.character(seu$tissue[seu$global_clusters==xx])[1])
  })
  return(factor(tis,levels=levels(seu$tissue)))
}

getTextColor <- function(hex) {
  ints <- strtoi(paste0('0x',c(substr(hex,2,3),substr(hex,4,5),substr(hex,6,7))))
  return({if (mean(ints)>100) { '#000000' } else {'#ffffff' }})
}

# ================= INIT DATA ========================

seu <- readRDS('./data/seurat_objects/rPCA_integrated.rds')

# ================= ANALYSIS =========================

# get centroids
cents <- sapply(levels(seu$global_clusters),function(cl){
  cells <- Cells(seu)[seu$global_clusters==cl]
  return(colMeans(Embeddings(seu,reduction='pca')[cells,]))
})
colnames(cents) <- levels(seu$global_clusters)
nj_tree <- nj(dist(t(cents)))
nj_tree$tissue <- sapply(nj_tree$tip.label,function(cl){ unname(seu$tissue[seu$global_clusters==cl][1]) })
plot.phylo(nj_tree)


# ================= PLOTTING =========================

p1 <- ggtree(nj_tree,layout='daylight',branch.length='none') 
p1$data$tissue <- getTissue(p1$data$label)
p1 <- p1 + geom_label(aes(label=label,fill=tissue),color='transparent') + 
  geom_text(aes(label=label,color=tissue)) + 
  scale_fill_manual(values=seu@misc$colors$tissue) +
  scale_color_manual(values=sapply(seu@misc$colors$tissue,getTextColor)) +
  coord_fixed() +
  guides(fill=guide_legend(title=NULL),color=guide_none())
ggsave('./figures/neighbor_joining_daylight.pdf',p1,width=8,height=8,units='in')

p2 <- ggtree(nj_tree) 
p2$data$tissue <- getTissue(p2$data$label)
p2 <- p2 + geom_label(aes(label=label,fill=tissue),color='transparent') + 
  geom_text(aes(label=label,color=tissue)) + 
  scale_fill_manual(values=seu@misc$colors$tissue) +
  scale_color_manual(values=sapply(seu@misc$colors$tissue,getTextColor)) +
  guides(fill=guide_legend(title=NULL),color=guide_none()) +
  theme(plot.title=NULL) +
  scale_y_continuous(expand=expansion(mult=c(.02,.02))) +
  scale_x_continuous(expand=expansion(mult=c(0,.05)))
ggsave('./figures/neighbor_joining_standard.pdf',p2,width=6,height=9,units='in')


