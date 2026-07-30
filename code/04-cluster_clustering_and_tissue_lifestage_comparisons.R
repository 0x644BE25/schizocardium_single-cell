######################################################
# CLUSTER CLUSTERING AND TISSUE LIFESTAGE COMPARISONS
#
# GOAL: Calculate cluster centroids (55 PC space) and
# perform hierarchical clustering.
#
# Look at endoderm+neural and muscle+mesoderm with
# cells separated by lifestage, and use inter-group
# PC centroids as well as cell-to-group PC distances
# to assess transcriptional similarity.
######################################################

setwd('/n/projects/cb2350/schizocardium_sc/')

# ================= IMPORTS ==========================

library(Seurat)
library(ggplot2)
library(patchwork)

# ================= PARAMS ===========================


# ================= METHODS ==========================

blendColors <- function(A,B,ratio) {
  ratio <- min(1,max(ratio,0))
  if (!startsWith(A,'#')) { A <- gplots::col2hex(A) }
  if (!startsWith(B,'#')) { B <- gplots::col2hex(B) }
  # ratio = b/total
  all <- NULL
  for (col in A) {
    #print(col)
    res <- '#'
    for (i in 1:3) {
      a <- strtoi(paste0('0x',substr(col,(2*i),2*i+1)))
      b <- strtoi(paste0('0x',substr(B,(2*i),2*i+1)))
      ab <- min(255,ceiling(((1-ratio)*a)+(ratio*b)))
      ab <- as.character(as.hexmode(ab))
      if (nchar(ab)==1) { ab <- paste0('0',ab)}
      res <- c(res,ab)
    }
    res <- paste0(res,collapse='')
    all <- c(all,res)
    #print(res)
  }
  return(all)
}

getTextColor <- function(A) {
  if (!startsWith(A,'#')) { A <- gplots::col2hex(A) }
  ints <- strtoi(paste0('0x',c(substr(A,2,3),substr(A,4,5),substr(A,6,7))))
  return({if (mean(ints)>128) { '#000000' } else {'#ffffff' }})
}

euc <- function(a,b) {
  return(sqrt(sum((a-b)^2)))
}

getCells <- function(seu,tissues=NA,clusters=NA,otherNotAll=TRUE) {
  cells <- list()
  
  # OUTGROUP
  if (otherNotAll) {
    cells[['larval_Other']] <- Cells(seu)[!seu$tissue %in% tissues & seu$timepoint %in% c('earlylarva','latelarva')]
    cells[['juvenile_Other']] <- Cells(seu)[!seu$tissue %in% tissues & seu$timepoint %in% c('earlyjuvenile','latejuvenile')]
  } else {
    cells[['larval_All']] <- Cells(seu)[seu$timepoint %in% c('earlylarva','latelarva')]
    cells[['juvenile_All']] <- Cells(seu)[seu$timepoint %in% c('earlyjuvenile','latejuvenile')]
  }
  
  # TISSUES
  for (tissue in setdiff(tissues,NA)) {
    cells[[paste0('larval_',tissue)]] <- Cells(seu)[seu$tissue==tissue & seu$timepoint %in% c('earlylarva','latelarva')]
    cells[[paste0('juvenile_',tissue)]] <- Cells(seu)[seu$tissue==tissue & seu$timepoint %in% c('earlyjuvenile','latejuvenile')]
  }
  
  # CLUSTERS
  for (cluster in setdiff(clusters,NA)) {
    cells[[paste0('cluster_',cluster)]] <- Cells(seu)[seu$global_clusters==cluster]
  }
  
  return(cells)
}

scaleToMaxIntercluster <- function(seu,curres) {
  maxdist <- max(dist(t(sapply(levels(seu$global_clusters),\(cl){
    colMeans(seu[['pca']]@cell.embeddings[Cells(seu)[seu$global_clusters==cl],])
  }))))
  
  curres$data$hcl$height <- 100*curres$data$hcl$height/maxdist
  curres$data$df_hm$value <- 100*curres$data$df_hm$value/maxdist
  curres$data$df_bp$dist <- 100*curres$data$df_bp$dist/maxdist
  
  return(curres)
}

larvJuvComparison <- function(seu,tissues=NA,clusters=NA,cols=NA,othercol='#666666',otherNotAll=TRUE,doScaling=TRUE) {
  library(ggplot2)
  library(gplots)
  library(ggdendro)
  library(ggh4x)
  library(patchwork)
  
  res <- list(data=list(PC_centroids=NA,hcl=NA,df_hm=NA,df_bp=NA),
              plots=list(dendrogram=NA,heatmap=NA,blank=NA,boxplot=NA))
  
  cells <- getCells(seu,tissues,clusters,otherNotAll=otherNotAll)
  
  # DISTANCES ====================
  emb <- Embeddings(seu,reduction='pca')
  PC_centroids <- t(sapply(cells,function(x){ colMeans(emb[x,]) }))
  res$data$PC_centroids <- PC_centroids
  hcl <- hclust(dist(PC_centroids))
  res$data$hcl <- hcl
  
  # FOR HEATMAP
  df_hm <- data.frame(as.matrix(dist(PC_centroids)))
  df_hm[df_hm==0] <- NA
  df_hm$y <- factor(rownames(df_hm),levels=rev(hcl$labels[hcl$order]))
  df_hm <- reshape2::melt(df_hm)
  df_hm$x <- factor(as.character(df_hm$variable),levels=rev(levels(df_hm$y)))
  res$data$df_hm <- df_hm
  
  # FOR BOXPLOTS
  df_bp <- do.call(rbind,lapply(names(cells),\(n){
    return(do.call(rbind,lapply(names(cells),\(n2){
      dists <- apply(emb[cells[[n]],],1,\(x){ euc(x,PC_centroids[n2,]) })
      return(data.frame(dist=unname(dists),cells=n,centroid=n2))
    })))
  }))
  df_bp$cells <- factor(as.character(df_bp$cells),levels=levels(df_hm$x))
  df_bp$centroid <- factor(df_bp$centroid,levels=levels(df_hm$x))
  res$data$df_bp <- df_bp
  
  roundTo <- 1
  if (doScaling) { res <- scaleToMaxIntercluster(seu,res); roundTo <- 0 }
  
  res$data$df_bp$x <- 0
  qs <- NULL
  for (cells in levels(res$data$df_bp$cells)) {
    for (cent in levels(res$data$df_bp$cent)) {
      res$data$df_bp[res$data$df_bp$cells==cells & res$data$df_bp$centroid==cent,'x'] <- rank(res$data$df_bp[res$data$df_bp$cells==cells & res$data$df_bp$centroid==cent,'dist'])/sum(res$data$df_bp$cells==cells & res$data$df_bp$centroid==cent)
      qs <- c(qs,quantile(res$data$df_bp[res$data$df_bp$cells==cells & res$data$df_bp$centroid==cent,'dist'],.95))
    }
  }
  
  lvs <- levels(res$data$df_bp$cells)
  #res$data$df2 <- do.call(rbind,lapply(lvs,\(cells){
  #  cell_dists <- res$data$df_bp[res$data$df_bp$cells==cells & res$data$df_bp$centroid==cells,'dist']
  #  cell_mean <- mean(cell_dists)
  #  do.call(rbind,lapply(lvs,\(cent){
  #    cent_dists <- res$data$df_bp[res$data$df_bp$cells==cells & res$data$df_bp$centroid==cent,'dist']
  #    ks <- ks.test(x=cent_dists,y=cell_dists)
  #    return(data.frame(cells=cells,centroid=cent,mean=tt$estimate,median=median(dists),stdev=sd(dists),conf.low=cell_mean+tt$conf.int[1],conf.high=cell_mean+tt$conf.int[2]))
  #  }))
  #}))
  
  # PLOTTING =====================
  
  if (is.na(cols[1])) {   cols <- c(seu@misc$colors$tissue) }
  if (!'Other' %in% names(cols)) { cols['Other'] <- othercol }
  if (!'All' %in% names(cols)) { cols['All'] <- othercol }
  
  cols <- c(setNames(cols,paste0('juvenile_',names(cols))),
            setNames(sapply(cols,\(x){ blendColors(x,'white',.5) }),paste0('larval_',names(cols))),
            setNames(rep('white',length(clusters)),paste0('cluster_',clusters)))
  
  txtcols <- sapply(cols,getTextColor)
  
  strip <- strip_themed(background_x=lapply(levels(df_hm$x),function(x){ element_rect(fill=cols[x]) }),
                        text_x=lapply(levels(df_hm$x),function(x){ element_text(color=txtcols[x]) }),
                        background_y=lapply(levels(df_hm$y),function(x){ element_rect(fill=cols[x]) }),
                        text_y=lapply(levels(df_hm$y),function(x){ element_text(color=txtcols[x]) }))
  
  stripLabel <- function(x){ sub('_','\n',x) }
  noLabel <- function(x){ return(NULL) }
  
  q <- max(unname(sapply(levels(res$data$df_bp$cells),\(cells){
    max(unname(sapply(levels(res$data$df_bp$centroid),\(cent){
      quantile(res$data$df_bp[res$data$df_bp$cells==cells & res$data$df_bp$centroid==cent,'dist'],.99)
    })))
  })))
  
  eb <- element_blank()
  ebtxt <- element_text(size=0,margin=margin(0,0,0,0,"cm"))
  
  res$data$df_bp$nest <- factor(gsub('_','\n',as.character(res$data$df_bp$cells)),levels=gsub('_','\n',levels(res$data$df_bp$cells)))
  strip2 <- strip_nested(background_x=c(lapply(levels(res$data$df_bp$cells),function(x){ element_rect(fill=cols[x]) }),rep(lapply(levels(res$data$df_bp$centroid),function(x){ eb }),length(levels(res$data$df_bp$centroid)))),
                         text_x=c(lapply(levels(res$data$df_bp$cells),function(x){ element_text(color=txtcols[x]) }),rep(lapply(levels(res$data$df_bp$centroid),function(x){ ebtxt }),length(levels(res$data$df_bp$centroid)))),
                         background_y=c(lapply(levels(res$data$df_bp$cells),function(x){ element_rect(fill=cols[x]) }),rep(lapply(levels(res$data$df_bp$centroid),function(x){ eb }),length(levels(res$data$df_bp$centroid)))),
                         text_y=c(lapply(levels(res$data$df_bp$cells),function(x){ element_text(color=txtcols[x]) }),rep(lapply(levels(res$data$df_bp$centroid),function(x){ ebtxt }),length(levels(res$data$df_bp$centroid)))))
  strip3 <- strip_nested(background_x=c(lapply(levels(res$data$df_bp$cells),function(x){ element_rect(fill=cols[x]) }),rep(lapply(levels(res$data$df_bp$centroid),function(x){ eb }),length(levels(res$data$df_bp$centroid)))),
                         text_x=c(lapply(levels(res$data$df_bp$cells),function(x){ element_text(color=txtcols[x]) }),rep(lapply(levels(res$data$df_bp$centroid),function(x){ ebtxt }),length(levels(res$data$df_bp$centroid)))),
                         background_y=c(lapply(levels(res$data$df_bp$cells),function(x){ element_rect(fill=cols[x]) }),rep(lapply(levels(res$data$df_bp$centroid),function(x){ eb }),length(levels(res$data$df_bp$centroid)))),
                         text_y=c(lapply(levels(res$data$df_bp$cells),function(x){ element_text(color=txtcols[x]) }),rep(lapply(levels(res$data$df_bp$centroid),function(x){ ebtxt }),length(levels(res$data$df_bp$centroid)))))
  
  # DENDROGRAM
  res$plots$dendrogram <- ggdendrogram(res$data$hcl,rotate=TRUE,size=0,labels=FALSE,leaf_labels=FALSE) + 
    scale_y_reverse() +
    theme(axis.text.y=element_blank(),
          axis.text.x=element_blank(),
          plot.margin=unit(c(.25,0,0,.25),'in'))
  
  # HEATMAP
  res$plots$heatmap <- ggplot(res$data$df_hm,aes(x=1,y=1,fill=value,label=round(value,roundTo),color=value>.7*max(value,na.rm=TRUE))) + 
    geom_tile(color='transparent',alpha=.8) + 
    facet_grid2(y~x,scales='free',space='free',switch='both',strip=strip,labeller=labeller(y=stripLabel,x=stripLabel)) +
    scale_fill_gradient(low='#dddddd',high='black',na.value='white') +
    geom_text(size=6) +
    scale_color_manual(values=c('black','white')) +
    theme_bw() +
    theme(panel.background=element_rect(color='white',fill='white'),
          panel.grid=element_blank(),axis.title=element_blank(),
          panel.spacing=unit(0,'in'),
          strip.text.y.left=element_text(angle=0),
          strip.text.x=element_text(angle=0),
          legend.position='none',
          plot.margin=unit(c(.25,.25,0,0),'in')) +
    scale_x_continuous(expand=c(0,0),breaks=NULL) +
    scale_y_continuous(expand=c(0,0),breaks=NULL)
  
  # SPACER
  res$plots$blank <- ggplot(data.frame()) + geom_blank() + theme_void() 
  
  # BOXPLOTS
  res$plots$boxplot <- ggplot(res$data$df_bp,aes(y=dist,fill=centroid,x=centroid)) + 
    geom_boxplot(outliers=FALSE,notch=TRUE,notchwidth=.2,linewidth=.25) +
    scale_fill_manual(values=cols) +
    facet_grid2(.~cells,switch='both',strip=strip,labeller=labeller(cells=stripLabel)) +
    theme_bw() +
    theme(axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          legend.position='none',
          axis.title=element_blank(),
          panel.grid.major.x=element_blank(),
          panel.grid.minor.y=element_blank(),
          panel.spacing=unit(0,'in'),
          plot.margin=unit(c(0,.25,.25,0),'in'))
  res$plots$boxplot <- patchwork::free(res$plots$boxplot,type='label',side='l')
  
  # CDF plots
  res$plots$cdf <- ggplot(res$data$df_bp,aes(x=x,y=dist,color=centroid,group=centroid)) + 
    geom_line(linewidth=.5) + 
    scale_color_manual(values=cols) +
    facet_grid2(.~cells,switch='both',strip=strip,labeller=labeller(cells=stripLabel)) +
    scale_y_continuous(trans='log2',limits=c(NA,max(qs)),name='distance to centroid') +
    scale_x_continuous(name=NULL) +
    theme_bw() + 
    theme(axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          legend.position='none',
          panel.grid.major.x=element_blank(),
          panel.grid.minor.x=element_blank(),
          panel.grid.minor.y=element_blank(),
          panel.spacing=unit(0,'in'),
          plot.margin=unit(c(0,.25,.25,0),'in'))
  res$plots$cdf <- patchwork::free(res$plots$cdf,type='label',side='l')
  
  # denstiy
  lvs <- factor(levels(res$data$df_bp$cells),levels=levels(res$data$df_bp$cells))
  nests <- factor(levels(res$data$df_bp$nest),levels=levels(res$data$df_bp$nest))
  df2 <- data.frame(nest=nests,centroid=lvs[1],x=0,xend=0)
  df2$y <- min(res$data$df_bp$dist); df2$yend <- q
  #df2h <- data.frame(nest=lvs[length(lvs)],
  #                   centroid=unique(res$data$df_bp$cells)[length(levels(res$data$df_bp$cells))])
  
  res$plots$dens <- ggplot(NULL) + 
    geom_density(data=res$data$df_bp,aes(y=dist,color=centroid,fill=centroid)) +
    facet_nested(.~nest+centroid,scales='free',strip=strip3,switch='both',labeller=labeller(nest=stripLabel)) +
    geom_vline(data=df2[1:6,],aes(xintercept=x),color='black',linewidth=1) +
    geom_hline(yintercept=0,color='black',linewidth=1) +
    geom_hline(yintercept=q,color='black',linewidth=1) +
    scale_fill_manual(values=cols) +
    scale_color_manual(values=cols) +
    theme_bw() +
    theme(panel.spacing=unit(0,'in'),
          panel.background=element_rect(color='#CCCCCC'),
          panel.border=element_blank(),
          panel.grid.major.x=element_blank(),
          panel.grid.minor.x=element_blank(),
          panel.grid.major.y=element_line(color='#CCCCCC',linetype='dashed'),
          axis.ticks.x=element_blank(),
          axis.text.x=element_blank(),
          axis.title.x=element_blank(),
          legend.position='none') +
    scale_y_continuous(trans='log2',limits=c(NA,q),expand=c(0,0),name='distance to centroid') +
    scale_x_continuous(limits=c(0,NA),expand=expansion(mult=c(0,.2)))
  res$plots$dens <- patchwork::free(res$plots$dens,type='label',side='l')
  return(res)
}

# ================= INIT DATA ========================

seu <- readRDS('./data/seurat_objects/rPCA_integrated.rds')
seu <- subset(seu,cells=Cells(seu)[seu$timepoint!='juvenileproboscis'])
seu$timepoint <- droplevels(seu$timepoint)

# ================= CLUSTER CLUSTERING ===============

emb <- Embeddings(seu,reduction='pca')
cents <- sapply(levels(seu$global_clusters),function(cl){
  cells <- Cells(seu)[seu$global_clusters==cl]
  colMeans(emb[cells,])
})
x <- as.matrix(dist(t(cents)))
colnames(x) <- rownames(x) <- paste0('cluster_',rownames(x))
write.csv('./data/inter_cluster_distances.csv',row.names=1)

cl <- hclust(dist(t(cents)))
df <- data.frame(cluster=cl$labels,order=cl$order,height=c(0,cl$height))
write.csv('./data/cluster_clustering_results.csv',row.names=FALSE)

# ================= LIFETSAVE VS TISSUE =========================

# DATA PREP
seu <- subset(seu,cells=Cells(seu)[seu$timepoint %in% c('earlylarva','latelarva','earlyjuvenile','latejuvenile')])
seu$timepoint <- droplevels(seu$timepoint)

# ENDODERM + NEURAL
comp <- larvJuvComparison(seu,tissues=c('Endoderm','Neural'),doScaling=TRUE)
p <- patchwork::wrap_plots(comp$plots[c('dendrogram','heatmap','blank','boxplot')],ncol=2) + plot_layout(widths=c(1,5),heights=c(5,2))
ggsave('./figures/lifestage_vs_endoderm_neural_boxplot.pdf',p,width=6,height=7,units='in',dpi=100)
p <- patchwork::wrap_plots(comp$plots[c('dendrogram','heatmap','blank','cdf')],ncol=2) + plot_layout(widths=c(1,5),heights=c(5,2))
ggsave('./figures/lifestage_vs_endoderm_neural_cdf.pdf',p,width=6,height=7,units='in',dpi=100)
p <- patchwork::wrap_plots(comp$plots[c('dendrogram','heatmap','blank','dens')],ncol=2) + plot_layout(widths=c(1,5),heights=c(5,2.25))
ggsave('./figures/lifestage_vs_endoderm_neural_density.pdf',p,width=6,height=7.25,units='in',dpi=100)

# MUSCLE + MESODERM
comp <- larvJuvComparison(seu,tissues=c('Muscle','Mesoderm'),doScaling=TRUE)
p <- patchwork::wrap_plots(comp$plots[c('dendrogram','heatmap','blank','boxplot')],ncol=2) + plot_layout(widths=c(1,5),heights=c(5,2))
ggsave('./figures/lifestage_vs_muscle_mesoderm_boxplot.pdf',p,width=6,height=7,units='in',dpi=100)
p <- patchwork::wrap_plots(comp$plots[c('dendrogram','heatmap','blank','cdf')],ncol=2) + plot_layout(widths=c(1,5),heights=c(5,2))
ggsave('./figures/lifestage_vs_muscle_mesoderm_cdf.pdf',p,width=6,height=7,units='in',dpi=100)
p <- patchwork::wrap_plots(comp$plots[c('dendrogram','heatmap','blank','dens')],ncol=2) + plot_layout(widths=c(1,5),heights=c(5,2.25))
ggsave('./figures/lifestage_vs_muscle_mesoderm_density.pdf',p,width=6,height=7.25,units='in',dpi=100)
