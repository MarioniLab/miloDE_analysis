
add_logcounts = function(sce, batch = "sample"){
  require(scran)
  require(scuttle)
  require(batchelor)
  clusters = quickCluster(sce, method="igraph", use.ranks=TRUE, d=50, min.mean=.1)
  sce = computeSumFactors(sce, clusters=clusters)

  meta = as.data.frame(colData(sce))
  batchFactor = factor(meta[, colnames(meta) == batch])
  sce = lapply(unique(batchFactor) , function(current.batch){
    idx = which(batchFactor == current.batch)
    current.sce = sce[,idx]
    clusters = quickCluster(current.sce, method="igraph", use.ranks=TRUE, d=50, min.mean=.1)
    current.sce = computeSumFactors(current.sce, clusters=clusters)
    return(current.sce)
  })
  sce = do.call(multiBatchNorm , sce )
  sce = do.call(cbind, sce)
  return(sce)
}

add_logcounts_simple = function(sce){
  require(scuttle)
  sce = logNormCounts(sce)
  return(sce)
}



add_batch_corrected_pca = function(sce, genes , reduced.dim_name = "pca.corrected" , batch = "sample" , d = 30, bpparam = NULL , correct = T){
  require(batchelor)
  #sce_ref = retain_informative_genes(sce)
  #hvgs = rownames(sce_ref)
  meta = as.data.frame(colData(sce))
  batchFactor = factor(meta[, colnames(meta) == batch])
  if (is.null(bpparam)){
    mbpca = multiBatchPCA(sce[genes , ], batch = batchFactor, d = d)
  }
  else {
    mbpca = multiBatchPCA(sce[genes , ], batch = batchFactor, d = d , BPPARAM = bpparam)
  }
  if (correct){
    out = do.call(reducedMNN, mbpca)
    joint.pca = out$corrected
  }
  else {
    joint.pca = do.call(rbind , mbpca)
  }
  joint.pca = joint.pca[order(rownames(joint.pca)) , ]
  sce = sce[, order(colnames(sce))]
  reducedDim(sce , reduced.dim_name) = joint.pca
  return(sce)
}



add_azimuth_pca = function(sce , ref_samples , query_samples, nPC = 30, reducedDim.name ,
                           bpparam){
  require(Seurat)
  require(SeuratObject)
  require(SingleCellExperiment)
  sce_seurat <- CreateSeuratObject(counts = counts(sce))
  sce_seurat = AddMetaData(sce_seurat, as.data.frame(colData(sce)), col.name = NULL)
  sce_seurat.list <- SplitObject(sce_seurat, split.by = "sample")
  # normalise and hvgs
  for (i in 1:length(sce_seurat.list)) {
    sce_seurat.list[[i]] <- NormalizeData(sce_seurat.list[[i]], verbose = FALSE)

    sce_seurat.list[[i]] <- FindVariableFeatures(sce_seurat.list[[i]], selection.method = "vst", nfeatures = 3000,
                                                 verbose = FALSE)
  }
  sce_reference.list = sce_seurat.list[ref_samples]
  # integrate reference
  ref_anchors <- FindIntegrationAnchors(object.list = sce_reference.list, anchor.features = 3000, dims = 1:nPC)
  sce_reference.integrated <- IntegrateData(anchorset = ref_anchors, dims = 1:nPC)
  DefaultAssay(sce_reference.integrated) <- "integrated"
  sce_reference.integrated <- ScaleData(sce_reference.integrated, verbose = FALSE)
  sce_reference.integrated <- RunPCA(sce_reference.integrated, npcs = nPC, verbose = FALSE)

  # map query
  pca_proj_query = bplapply(query_samples , function(current.sample){
    current.sce = sce_seurat.list[[which(names(sce_seurat.list) == current.sample)]]
    query_anchors <- FindTransferAnchors(reference = sce_reference.integrated, query = current.sce,
                                         dims = 1:nPC, reference.reduction = "pca")
    predictions <- TransferData(anchorset = query_anchors, refdata = t(Embeddings(sce_reference.integrated[['pca']])), dims = 1:nPC)
    predictions = predictions[1:nrow(predictions), 1:ncol(predictions)]
    current.sce[["pca"]] <- CreateDimReducObject(embeddings = t(as.matrix(predictions)), key = "PC_", assay = DefaultAssay(current.sce))
    current.sce = as.SingleCellExperiment(current.sce)
    out = reducedDim(current.sce , "PCA")
    colnames(out) = paste0("PC_" , c(1:nPC))
    return(out)
  } , BPPARAM = bpparam)

  # combine
  pca_proj_query = do.call(rbind , pca_proj_query)
  sce_reference.integrated = as.SingleCellExperiment(sce_reference.integrated)
  pca_ref = reducedDim(sce_reference.integrated , "PCA")
  colnames(pca_ref) = paste0("PC_" , c(1:nPC))
  pca_joint = rbind(pca_proj_query , pca_ref)
  pca_joint = pca_joint[order(rownames(pca_joint)) , ]
  sce = sce[, order(colnames(sce))]
  reducedDim(sce , reducedDim.name) = pca_joint

  #umaps = as.data.frame(uwot::umap(reducedDim(sce , reducedDim.name) , min_dist = 0.7))
  #rownames(umaps) = colnames(sce)
  #colnames(umaps) = c("x" , "y")
  #reducedDim(sce , umap.name) = umaps

  return(sce)
}



retain_informative_genes = function(sce, n = NULL, var.thresh = 0){
  require(scran)
  dec.sce = modelGeneVar(sce)
  hvg.genes = getTopHVGs(dec.sce, var.threshold = var.thresh)
  if (!is.null(n)){
    if (n < length(hvg.genes)){
      hvg.genes = getTopHVGs(dec.sce, n = n)
    }
  }
  sce = sce[hvg.genes, ]
  return(sce)
}


#gene_df is the cellranger gene table
getHVGs = function(sce, var.thresh = 0, n = NULL){
  dec.sce = modelGeneVar(sce)
  hvg.genes = getTopHVGs(dec.sce, var.threshold = var.thresh)
  if (!is.null(n)){
    if (n < length(hvg.genes)){
      hvg.genes = getTopHVGs(dec.sce, n = n)
    }
  }
  return(hvg.genes)
}


getmode <- function(v, dist) {
  tab = table(v)
  #if tie, break to shortest distance
  if(sum(tab == max(tab)) > 1){
    tied = names(tab)[tab == max(tab)]
    sub = dist[v %in% tied]
    names(sub) = v[v %in% tied]
    return(names(sub)[which.min(sub)])
  } else {
    return(names(tab)[which.max(tab)])
  }
}


cols_order = c("#3b7c70", "#ce9642")
names(cols_order) = c(1,2)



