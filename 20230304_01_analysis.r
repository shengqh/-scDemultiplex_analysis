source('/home/shengq2/program/scDemultiplex_analysis/common.r')

sample_map = list(
  "hto12"=list("HTO"="hto12_hto_mtx.rds", "RNA"="hto12_umi_mtx.rds"),
  "pbmc"=list("HTO"="pbmc_hto_mtx.rds", "RNA"="pbmc_umi_mtx.rds")
)

ignored_tags<-c("no_match","ambiguous","total_reads","bad_struct")

save_to_matrix<-function(counts, target_folder) {
  if(!dir.exists(target_folder)){
    dir.create(target_folder)
  }
  
  bar_file=paste0(target_folder, "/barcodes.tsv")
  writeLines(colnames(counts), bar_file)
  gzip(bar_file, overwrite=T)
  
  feature_file=paste0(target_folder, "/features.tsv")
  writeLines(rownames(counts), feature_file)
  gzip(feature_file, overwrite=T)
  
  matrix_file=paste0(target_folder, "/matrix.mtx")
  writeMM(counts, matrix_file)
  gzip(matrix_file, overwrite=T)
}

prepare_data<-function(root_dir, sample_map, cur_sample) {
  sample_folder=paste0(root_dir, cur_sample)
  if(!dir.exists(sample_folder)){
    dir.create(sample_folder)
  }
  setwd(sample_folder)

  hto_file=paste0(cur_sample, "_hto_mtx.rds")
  exp_file=paste0(cur_sample, "_umi_mtx.rds")

  if(!file.exists(hto_file)){
    smap = sample_map[[cur_sample]]
    hto_mtx_file=smap$HTO
    exp_mtx_file=smap$RNA
    hto=data.frame(fread(hto_mtx_file), row.names=1)

    #transpose the matrix
    if(nrow(hto) > 100){
      hto=data.frame(t(hto), check.names = F)
    }

    hto<-hto[!(rownames(hto) %in% ignored_tags), ]
    rownames(hto)<-gsub("-.+", "", rownames(hto))
    saveRDS(hto, hto_file)

    if(!is.null(exp_mtx_file)){    
      exp=data.frame(fread(exp_mtx_file), row.names=1)
      if(all(grepl("-",rownames(exp)))){
        exp=data.frame(t(exp), check.names = F)
      }
      saveRDS(exp, exp_file)
    }
  }else{
    hto=readRDS(hto_file)
    if(cur_sample == "hto12"){
      hto=t(hto)
    }
    if(file.exists(exp_file)){
      exp=readRDS(exp_file)
    }
  }

  if(file.exists(exp_file)){
    common_cells=colnames(hto)[colnames(hto) %in% colnames(exp)]
    print(paste0("common_cells = ", length(common_cells)))
    exp=exp[,common_cells]
    
    hto.exp <- CreateSeuratObject(counts = exp, min.features = 300)
    cells.valid<-colnames(hto.exp)
    hto<-hto[!(rownames(hto) %in% ignored_tags), cells.valid]
    
    DefaultAssay(hto.exp)<-"RNA"
    # Normalize RNA data with log normalization
    hto.exp <- NormalizeData(hto.exp)
    # Find and scale variable features
    hto.exp <- FindVariableFeatures(hto.exp, selection.method = "mean.var.plot")
    features = VariableFeatures(hto.exp)
    hto.exp <- ScaleData(hto.exp, features = features)    
    hto.exp <- RunPCA(hto.exp, features = features)
    hto.exp <- RunUMAP(hto.exp, reduction = "pca", dims=1:20)
    hto.exp <- FindNeighbors(hto.exp, reduction = "pca")
    hto.exp <- FindClusters(hto.exp, reduction = "pca", resolution=0.02) 
    
    hto.exp[["HTO"]] <- CreateAssayObject(counts = hto)
    saveRDS(hto.exp, paste0(cur_sample, ".combined.rds"))
  }
  counts=as.sparse(hto)
  
  save_to_matrix(counts=counts, target_folder="data")
  
  rds_file=paste0(cur_sample, ".counts.rds")
  saveRDS(counts, rds_file)
  
  obj <- scDemultiplex:::read_hto(rds_file)
  obj<-hto_umap(obj)
  saveRDS(obj, paste0(cur_sample, ".obj.rds"))
}

do_bff_raw<-function(root_dir, cur_sample) {
  library(cellhashR)

  sample_folder=paste0(root_dir, cur_sample)
  setwd(sample_folder)

  final_rds = paste0( "bff_raw/", cur_sample, ".bff_raw.rds")
  if(!file.exists(final_rds)){
    obj<-readRDS(paste0(cur_sample, ".obj.rds"))
    
    barcodeData <- obj$HTO@counts

    result_folder=paste0(sample_folder, "/bff_raw")
    if(!dir.exists(result_folder)){
      dir.create(result_folder)
    }
    setwd(result_folder)  

    tic(paste0("starting ", cur_sample, "...\n"))
    bff_calls <- GenerateCellHashingCalls(barcodeMatrix = barcodeData, 
                                              methods = c('bff_raw'),
                                              doTSNE = FALSE,
                                              doHeatmap = FALSE,
                                              bff_raw.min_average_reads = 0)  
    toc1=toc()
    
    rownames(bff_calls)<-bff_calls$cellbarcode
    bff_calls<-bff_calls[colnames(obj),]
    
    stopifnot(all(colnames(obj) == bff_calls$cellbarcode))
    
    obj$bff_raw <- bff_calls$bff_raw
    obj$bff_raw.global <- as.character(bff_calls$bff_raw)
    obj$bff_raw.global[!(obj$bff_raw.global %in% c("Negative", "Doublet"))] <- "Singlet"

    saveRDS(list("bff_raw"=toc1), paste0(cur_sample, ".bff_raw.tictoc.rds"))
    
    obj<-hto_plot(obj, paste0(cur_sample, ".bff_raw"), group.by="bff_raw")

    saveRDS(obj, paste0(cur_sample, ".bff_raw.rds"))
  }
}

do_bff_cluster<-function(root_dir, cur_sample) {
  library(cellhashR)

  sample_folder=paste0(root_dir, cur_sample)
  setwd(sample_folder)
  
  final_rds = paste0( "bff_cluster/", cur_sample, ".bff_cluster.rds")
  if(!file.exists(final_rds)){
    obj<-readRDS(paste0(cur_sample, ".obj.rds"))
    
    barcodeData <- obj$HTO@counts
    
    stopifnot(all(colnames(obj) == colnames(barcodeData)))
    
    result_folder=paste0(sample_folder, "/bff_cluster")
    if(!dir.exists(result_folder)){
      dir.create(result_folder)
    }
    setwd(result_folder)  
    
    tic(paste0("starting ", cur_sample, "...\n"))
    bff_calls <- GenerateCellHashingCalls(barcodeMatrix = barcodeData, 
                                              methods = c('bff_cluster'),
                                              doTSNE = FALSE,
                                              doHeatmap = FALSE,
                                              bff_cluster.min_average_reads = 0)  
    toc1=toc()
    
    rownames(bff_calls)<-bff_calls$cellbarcode
    bff_calls<-bff_calls[colnames(obj),]
    
    stopifnot(all(colnames(obj) == bff_calls$cellbarcode))
    
    obj$bff_cluster <- bff_calls$bff_cluster
    obj$bff_cluster.global <- as.character(bff_calls$bff_cluster)
    obj$bff_cluster.global[!(obj$bff_cluster.global %in% c("Negative", "Doublet"))] <- "Singlet"
    
    saveRDS(list("bff_cluster"=toc1), paste0(cur_sample, ".bff_cluster.tictoc.rds"))
    
    obj<-hto_plot(obj, paste0(cur_sample, ".bff_cluster"), group.by="bff_cluster")
    
    saveRDS(obj, paste0(cur_sample, ".bff_cluster.rds"))
  }
}

do_GMMDemux<-function(root_dir, sample_tags, cur_sample){
  sample_folder=paste0(root_dir, cur_sample)
  tags = sample_tags[[cur_sample]]
  setwd(sample_folder)

  final_rds = paste0( "GMM-demux/GMM_full.csv")
  if(!file.exists(final_rds)){
    cmd = paste0("/data/cqs/softwares/conda_py3_10/bin/GMM-demux data ", tags,  " -f GMM-demux -o GMM-demux")
    print(cmd)
    system(cmd)
  }
}

do_scDemultiplex<-function(root_dir, cur_sample, p.cut=0.001){
  sample_folder=paste0(root_dir, cur_sample)
  setwd(sample_folder)
  
  final_rds = paste0( "scDemultiplex/", cur_sample, ".scDemultiplex.rds")
  if(!file.exists(final_rds)){
    rds_file=paste0(cur_sample, ".obj.rds")
    obj=readRDS(rds_file)

    ntags = nrow(obj)

    result_folder=paste0(sample_folder, "/scDemultiplex")
    if(!dir.exists(result_folder)){
      dir.create(result_folder)
    }
    setwd(result_folder)

    final_file=paste0(cur_sample, ".scDemultiplex.rds")
    
    output_prefix<-paste0(cur_sample, ".HTO")

    tic(paste0("starting ", cur_sample, " cutoff ...\n"))
    cat("  scDemultiplex_cutoff ...\n")
    obj<-demulti_cutoff(obj, output_prefix, cutoff_startval = 0, mc.cores=ntags)
    toc1=toc()
    cat("  scDemultiplex_full ...\n")
    obj<-demulti_refine(obj, p.cut, refine_negative_doublet_only=FALSE, mc.cores=ntags)
    obj$scDemultiplex_full=obj$scDemultiplex
    obj$scDemultiplex_full.global=obj$scDemultiplex.global
    toc3=toc()
    cat("  scDemultiplex_negative_doublet_only ...\n")
    obj<-demulti_refine(obj, p.cut, refine_negative_doublet_only=TRUE, mc.cores=ntags)
    toc2=toc()

    saveRDS(list("cutoff"=toc1, "refine"=toc2, "full"=toc3), paste0(cur_sample, ".scDemultiplex.tictoc.rds"))

    obj<-hto_plot(obj, paste0(output_prefix, ".cutoff"), group.by="scDemultiplex_cutoff")
    obj<-hto_plot(obj, paste0(output_prefix, ".refine_p", p.cut), group.by="scDemultiplex")
    obj<-hto_plot(obj, paste0(output_prefix, ".full_p", p.cut), group.by="scDemultiplex_full")

    saveRDS(obj, final_file)
  }
}

do_seurat_HTODemux<-function(root_dir, cur_sample){
  sample_folder=paste0(root_dir, cur_sample)
  setwd(sample_folder)

  final_rds = paste0( "HTODemux/", cur_sample, ".HTODemux.rds")
  if(!file.exists(final_rds)){
    rds_file=paste0(cur_sample, ".obj.rds")
    obj=readRDS(rds_file)
    
    result_folder=paste0(sample_folder, "/HTODemux")
    if(!dir.exists(result_folder)){
      dir.create(result_folder)
    }
    setwd(result_folder)

    tic(paste0("starting ", cur_sample, "...\n"))
    obj <- HTODemux(obj, assay = "HTO", positive.quantile = 0.99)
    toc1=toc()

    obj$HTODemux <- obj$HTO_classification
    obj$HTODemux[which(obj$HTO_classification.global == "Doublet")] <- "Doublet"
    
    saveRDS(list("HTODemux"=toc1), paste0(cur_sample, ".HTODemux.tictoc.rds"))
    
    obj<-hto_plot(obj, paste0(cur_sample, ".HTODemux"), group.by="HTODemux")

    saveRDS(obj, paste0(cur_sample, ".HTODemux.rds"))
  }
}

# ----------------------------------------------------------------
# 3. the heuristic classifier of MULTI-seq - R
# https://github.com/chris-mcginnis-ucsf/MULTI-seq
# https://satijalab.org/seurat/reference/multiseqdemux

do_seurat_MULTIseqDemux<-function(root_dir, cur_sample){
  sample_folder=paste0(root_dir, cur_sample)
  setwd(sample_folder)
  
  final_rds = paste0( "MULTIseqDemux/", cur_sample, ".MULTIseqDemux.rds")
  if(!file.exists(final_rds)){
    rds_file=paste0(cur_sample, ".obj.rds")
    obj=readRDS(rds_file)
    
    result_folder=paste0(sample_folder, "/MULTIseqDemux")
    if(!dir.exists(result_folder)){
      dir.create(result_folder)
    }
    setwd(result_folder)
    
    tagnames<-rownames(obj)
    
    tic(paste0("starting ", cur_sample, "...\n"))
    obj <- MULTIseqDemux(obj, assay = "HTO")
    toc1=toc()
    
    obj$MULTIseqDemux <- as.character(obj$MULTI_classification)
    obj$MULTIseqDemux[which(!obj$MULTI_classification %in% c(tagnames, "Negative"))] <- "Doublet"
    
    saveRDS(list("MULTIseqDemux"=toc1), paste0(cur_sample, ".MULTIseqDemux.tictoc.rds"))
    
    obj<-hto_plot(obj, paste0(cur_sample, ".MULTIseqDemux"), group.by="MULTIseqDemux")
    
    saveRDS(obj, paste0(cur_sample, ".MULTIseqDemux.rds"))
  }
}

do_analysis<-function(root_dir, sample_map, sample_tags, cur_sample){
  obj_file = paste0(root_dir, "/", cur_sample, "/", cur_sample, ".obj.rds")
  if(!file.exists(obj_file)){
    cat("prepare_data ...\n")
    prepare_data(root_dir, sample_map, cur_sample)
  }

  bff_raw_file = paste0(root_dir, "/", cur_sample, "/bff_raw/", cur_sample, ".bff_raw.rds")
  if(!file.exists(bff_raw_file)){
    cat("bff_raw ...\n")
    do_bff_raw(root_dir, cur_sample)
  }

  bff_cluster_file = paste0(root_dir, "/", cur_sample, "/bff_cluster/", cur_sample, ".bff_cluster.rds")
  if(!file.exists(bff_cluster_file)){
    cat("bff_cluster ...\n")
    do_bff_cluster(root_dir, cur_sample)
  }

  gmm_file = paste0(root_dir, "/", cur_sample, "/GMM-demux/GMM_full.csv")
  if(!file.exists(gmm_file)){
    cat("GMM-demux ...\n")
    do_GMMDemux(root_dir, sample_tags, cur_sample)
  }

  scDemultiplex_file = paste0(root_dir, "/", cur_sample, "/scDemultiplex/", cur_sample, ".scDemultiplex.rds")
  if(!file.exists(scDemultiplex_file)){
    cat("scDemultiplex ...\n")
    do_scDemultiplex(root_dir, cur_sample)
  }

  htodemux_file = paste0(root_dir, "/", cur_sample, "/HTODemux/", cur_sample, ".HTODemux.rds")
  if(!file.exists(htodemux_file)){
    cat("HTODemux ...\n")
    do_seurat_HTODemux(root_dir, cur_sample)
  }

  MULTIseqDemux_file = paste0(root_dir, "/", cur_sample, "/MULTIseqDemux/", cur_sample, ".MULTIseqDemux.rds")
  if(!file.exists(MULTIseqDemux_file)){
    cat("MULTIseqDemux ...\n")
    do_seurat_MULTIseqDemux(root_dir, cur_sample)
  }

  cat("done.\n")
}

cur_sample = "hto12"
for(cur_sample in samples){
  do_analysis(root_dir, sample_map, sample_tags, cur_sample)
}
