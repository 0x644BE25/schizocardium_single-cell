######################################################
# NORMALIZE DIMREDUC CLUSTER
#
# GOAL: Prep, integrate, and perform dimensional
# reduction and clustering on all three batches of 
# schizocardium scRNA-seq expression matrices.
######################################################

# ================= IMPORTS ==========================

library('Seurat')
library('Matrix')  # mtx format

# ================= PARAMS ===========================

buildS3objects <- TRUE
rPCA <- TRUE

maxGB <- 40
nCores <- 16
options(future.globals.maxSize=(1024^3)*maxGB)
options(future.seed=TRUE) 

nPCs <- 55
nNeighbors <- 80

saveAs <- '~/schizocardium/data/seuratObjects/'
batches <- list(c('earlylarva','juvenileproboscis'),
                c('latelarva','duringmeta','earlyjuvenile','latejuvenile'),
                c('latelarva','duringmeta','earlyjuvenile','latejuvenile'))
allSamples <- c('earlylarva','latelarva','earlymeta','duringmeta','earlyjuvenile','latejuvenile','juvenileproboscis')
sampleLevels <- c('earlylarva_1',
                  'latelarva_2','latelarva_3',
                  'duringmeta_2','duringmeta_3',
                  'earlyjuvenile_2','earlyjuvenile_3',
                  'latejuvenile_2','latejuvenile_3',
                  'juvenileproboscis_1')
abbv <- list('earlylarva'='EL','latelarva'='LL','duringmeta'='DM','earlyjuvenile'='EJ','latejuvenile'='LJ','juvenileproboscis'='JP')
filePath <- '.data/expression_matrices/'
batch3suffix <- '_11062020'
scLevels <- paste(rep(abbv,each=100),0:99,sep='-')

# ================= BUILD SEURAT OBJECTS =============

if (buildS3objects) {
  for (batch in 1:length(batches)) {
    print(batch)
    samples <- batches[[batch]]
    suffix <- if (batch==3) { batch3suffix } else { '' }
    
    all <- NULL
    for (sample in samples){
      print(sample)
      mat <- readMM(paste0(filePath,sample,suffix,'_matrix.mtx'))
      feat <- read.csv(paste0(filePath,sample,suffix,'_features.tsv'),sep='\t',header=FALSE,stringsAsFactors=FALSE,col.names=c('SCA','Alias','Data'))
      barc <- read.csv(paste0(filePath,sample,suffix,'_barcodes.tsv'),sep='\t',header=FALSE,stringsAsFactors=FALSE)[,1]
      barc <- paste0(abbv[[sample]],'-',batch,'-',barc)
      rownames(mat) <- feat[,'Alias']
      colnames(mat) <- barc
      
      #curr[['RNA']]@meta.features$SCA <- feat$SCA
      if (is.null(all)) { 
        all <- mat
      } else {
        all <- cbind(all,mat)
      }
    }
    curr <- CreateSeuratObject(counts=all,project='schizo')
    curr$timepoint <- setNames(factor(setNames(names(abbv),abbv)[substr(Cells(curr),1,2)],levels=allSamples),Cells(curr))
    curr$batch <- factor(batch,levels=1:4)
    curr <- SCTransform(curr,variable.features.n=3000,return.only.var.genes=FALSE)
    saveRDS(curr,paste0('./data/seurat_objects/','batch_',batch,'.rds'))
  }
}

# ================= DO PCA ==================

if (rPCA) {
  seuratList <- list(readRDS(paste0('./data/seurat_objects/','batch_1.rds')),
                     readRDS(paste0('./data/seurat_objects/','batch_2.rds')),
                     readRDS(paste0('./data/seurat_objects/','batch_3.rds')))
  seuratList <- lapply(seuratList, function(x) {
    x@active.assay <- 'RNA'
    x <- NormalizeData(x)
    x <- FindVariableFeatures(x,selection.method="vst",nfeatures=3000)
  })
  
  allfeat <- unique(unlist(lapply(seuratList,function(x){rownames(x[['RNA']]$data)})))
  
  features <- SelectIntegrationFeatures(object.list=seuratList)
  seuratList <- lapply(seuratList,function(x) {
    x <- ScaleData(x,features=features,verbose=FALSE)
    x <- RunPCA(x,features=features,verbose=FALSE,npcs=nPCs)
  })
  
  anchors <- FindIntegrationAnchors(seuratList,anchor.features=features,reduction="rpca")
  integrated <- IntegrateData(anchorset=anchors,features.to.integrate=allfeat)
  saveRDS(integrated,paste0('./data/seurat_objects/','rPCA_integrated.rds'))
  DefaultAssay(integrated) <- 'integrated'
  integrated <- ScaleData(integrated)
  integrated <- RunPCA(integrated,npcs=55)
  integrated <- RunUMAP(integrated,dims=1:55,n.neighbors=80)
  integrated <- FindNeighbors(integrated,dims=1:55)
  integrated <- FindClusters(integrated)
  integrated$global_clusters <- factor(integrated$seurat_clusters)
  integrated$sample <- paste(integrated$timepoint,integrated$batch,sep='_')
  integrated$sample <- factor(integrated$sample,levels=sampleLevels)
  saveRDS(integrated,paste0('./data/seurat_objects/','rPCA_integrated.rds'))
}
