rm(list=ls()) 

source("/home/shengq2/program/scDemultiplex_analysis/common.r")

for(cur_sample in samples) {
  cat("Processing ", cur_sample, "...\n")
  setwd(paste0(root_dir, cur_sample))
  
  obj<-readRDS(paste0("scDemultiplex/", cur_sample, ".scDemultiplex.rds"))
  meta<-obj@meta.data
  meta$scDemultiplex = meta$scDemultiplex_full
  meta$scDemultiplex.global = meta$scDemultiplex_full.global
  
  HTODemux_meta<-readRDS(paste0("HTODemux/", cur_sample, ".HTODemux.rds"))@meta.data
  stopifnot(nrow(HTODemux_meta) == nrow(meta))
  HTODemux_meta<-HTODemux_meta[rownames(meta),]
  
  MULTIseqDemux_meta<-readRDS(paste0("MULTIseqDemux/", cur_sample, ".MULTIseqDemux.rds"))@meta.data
  stopifnot(nrow(MULTIseqDemux_meta) == nrow(meta))
  MULTIseqDemux_meta<-MULTIseqDemux_meta[rownames(meta),]
  
  dd3 <- read.csv("GMM-demux/GMM_full.csv", stringsAsFactors = F)
  dd3c <- read.csv("GMM-demux/GMM_full.config", header=FALSE, stringsAsFactors = F)
  gmm <- merge(dd3, dd3c, by.x = "Cluster_id", by.y = "V1", all.x = T, all.y = F)
  rownames(gmm)<-gmm$X
  gmm$GMM_demux <- gmm$V2
  gmm$GMM_demux <- gsub("\\s+","",gmm$GMM_demux)
  gmm$GMM_demux <- gsub("_", "-", gmm$GMM_demux)
  gmm$GMM_demux <- gsub("negative", "Negative", gmm$GMM_demux)
  gmm$GMM_demux[which(! gmm$GMM_demux %in% c("Negative", rownames(obj)))] <- "Doublet"
  stopifnot(nrow(gmm) == nrow(meta))
  gmm<-gmm[rownames(meta),]

  bff_raw_meta<-readRDS(paste0("bff_raw/", cur_sample, ".bff_raw.rds"))@meta.data
  bff_raw_meta$BFF_raw <- gsub("_", "-", bff_raw_meta$bff_raw)
  stopifnot(nrow(bff_raw_meta) == nrow(meta))
  bff_raw_meta<-bff_raw_meta[rownames(meta),]
  
  bff_cluster_meta<-readRDS(paste0("bff_cluster/", cur_sample, ".bff_cluster.rds"))@meta.data
  bff_cluster_meta$BFF_cluster <- gsub("_", "-", bff_cluster_meta$bff_cluster)
  stopifnot(nrow(bff_cluster_meta) == nrow(meta))
  bff_cluster_meta<-bff_cluster_meta[rownames(meta),]
  
  meta$HTODemux = factor(HTODemux_meta$HTODemux, levels=levels(obj$scDemultiplex))
  meta$MULTIseqDemux = factor(MULTIseqDemux_meta$MULTIseqDemux, levels=levels(obj$scDemultiplex))
  meta$GMM_Demux = factor(gmm$GMM_demux, levels=levels(obj$scDemultiplex))
  meta$BFF_raw = factor(bff_raw_meta$BFF_raw, levels=levels(obj$scDemultiplex))
  meta$BFF_cluster = factor(bff_cluster_meta$BFF_cluster, levels=levels(obj$scDemultiplex))
  
  write.csv(meta, file = paste0(cur_sample, ".results.csv"), row.names = T)
  
  obj@meta.data<-meta

  final_file = paste0(cur_sample, ".results_obj.rds")
  saveRDS(obj, final_file)
}
