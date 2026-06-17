######################################################
# NORMALIZATION, DIMENSIONAL REDUCTION, CLUSTERING
#
# GOAL: Create Seurat objects, prep and integrate all 
# three batches of schizocardium scRNA-seq expression 
# matrices.
######################################################

# ================= IMPORTS ==========================

library(Seurat)

# ================= PARAMS ===========================

buildS3objects <- TRUE
rPCAintegrate <- TRUE
doBasicAnalysis <- TRUE

maxGB <- 64
options(future.globals.maxSize=(1024^3)*maxGB)
options(future.seed=TRUE) 

nPCs <- 55
nNeighbors <- 80

readDir <- './data/expression_matrices/'
writeDir <- './data/seurat_objects/'

batches <- list(c('earlylarva','juvenileproboscis'),
                c('latelarva','duringmeta','earlyjuvenile','latejuvenile'),
                c('latelarva','duringmeta','earlyjuvenile','latejuvenile'))
batch3suffix <- '_11062020'
abbv <- list('earlylarva'='EL','latelarva'='LL',
             'duringmeta'='DM',
             'earlyjuvenile'='EJ','latejuvenile'='LJ',
             'juvenileproboscis'='JP')

# ================= BUILD SEURAT OBJECTS =============

if (buildS3objects) {
  if (!dir.exists(writeDir)) { dir.create(writeDir) }
  
  seuList <- list()
  for (batch in 1:length(batches)) {
    print(batch)
    samples <- batches[[batch]]
    
    mat <- NULL
    for (sample in samples){
      print(sample)
      curr_mat <- ReadMtx(paste0(readDir,'batch',batch,'_',sample,'_matrix.mtx'),
                          cells=paste0(readDir,'batch',batch,'_',sample,'_barcodes.tsv'),cell.column=1,cell.sep='\t',
                          features=paste0(readDir,'batch',batch,'_',sample,'_features.tsv'),feature.column=2,feature.sep='\t')
      colnames(curr_mat) <- paste(paste0(abbv[[sample]],'-',batch,'-',colnames(curr_mat)))
      mat <- cbind(mat,curr_mat)
    }
    seu <- CreateSeuratObject(counts=mat,project='schizo_sc')
    seu <- SCTransform(seu,ncells=length(Cells(seu)),variable.features.n=3000,vst.flavor='v1',return.only.var.genes=FALSE)
    saveRDS(seu,paste0('./data/seurat_objects/batch_',batch,'.rds'),compress=FALSE)
  }
}
# ================= RPCA INTEGRATION =================

if (rPCAintegrate) {
  seuList <- lapply(1:3,\(i){ readRDS(paste0('./data/seurat_objects/batch_',i,'.rds')) })
  seuList <- lapply(seuList,function(x) {
    x@active.assay <- 'RNA'
    x <- NormalizeData(x)
    x <- FindVariableFeatures(x,selection.method="vst",nfeatures=3000)
  })
  
  allfeat <- unique(unlist(lapply(seuList,function(x){ rownames(x[['RNA']]) })))
  
  features <- SelectIntegrationFeatures(object.list=seuList)
  seuList <- lapply(seuList,function(x) {
    x <- ScaleData(x,features=allfeat,verbose=FALSE)
    x <- RunPCA(x,features=features,verbose=FALSE,npcs=nPCs)
  })
  
  anchors <- FindIntegrationAnchors(seuList,anchor.features=features,reduction="rpca")
  rm(seuList)
  
  integrated <- IntegrateData(anchors,features.to.integrate=allfeat)
}

# ================= BASIC ANALYSIS ===================

if (doBasicAnalysis) {
  DefaultAssay(integrated) <- 'integrated'
  integrated <- ScaleData(integrated)
  integrated <- RunPCA(integrated,npcs=nPCs)
  integrated <- RunUMAP(integrated,dims=1:nPCs,n.neighbors=nNeighbors)
  integrated <- FindNeighbors(integrated,dims=1:nPCs)
  integrated <- FindClusters(integrated)
  
  saveRDS(integrated,paste0(writeDir,'rPCA_integrated.rds'),compress=FALSE)
}
