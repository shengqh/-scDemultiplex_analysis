rm(list=ls()) 

source("/home/shengq2/program/scDemultiplex_analysis/common.r")

cur_sample=samples[3]
for(cur_sample in samples) {
  cat("Processing ", cur_sample, "...\n")
  setwd(paste0(root_dir, cur_sample))
  
  scDemultiplex.p.cut = scDemultiplex.p.cuts[1]
  rds_file = paste0(get_scDemultiplex_folder(root_dir, cur_sample, scDemultiplex.p.cut), "/", cur_sample, ".scDemultiplex.rds")
  obj<-readRDS(rds_file)
  meta<-obj@meta.data
  meta$scDemultiplex = as.character(meta$scDemultiplex_full)
  meta$scDemultiplex.global = meta$scDemultiplex_full.global
  
  detail_file = paste0(get_scDemultiplex_folder(root_dir, cur_sample, scDemultiplex.p.cut), "/", cur_sample, ".HTO.iteration.detail.csv")
  if(file.exists(detail_file)){
    details = read.csv(detail_file, header=T, row.names=1)
    for(i in 1:ncol(details)){
      details[,i] <- as.character(details[,i])
    }
    stopifnot(all.equal(rownames(details), rownames(meta)))
    meta<-cbind(meta, details)
  }
  
  if(length(scDemultiplex.p.cuts) > 1){
    for(idx in c(2:length(scDemultiplex.p.cuts))){
      scDemultiplex.p.cut = scDemultiplex.p.cuts[idx]
      scDemultiplex_meta<-readRDS(paste0(get_scDemultiplex_folder(root_dir, cur_sample, scDemultiplex.p.cut), "/", cur_sample, ".scDemultiplex.rds"))@meta.data
      stopifnot(nrow(scDemultiplex_meta) == nrow(meta))
      scDemultiplex_meta<-scDemultiplex_meta[rownames(meta),]
      s_name = names(scDemultiplex.p.cuts)[idx]
      meta[[s_name]] = as.character(scDemultiplex_meta$scDemultiplex_full)
    }
  }
  
  HTODemux_meta<-readRDS(paste0("HTODemux/", cur_sample, ".HTODemux.rds"))@meta.data
  stopifnot(nrow(HTODemux_meta) == nrow(meta))
  HTODemux_meta<-HTODemux_meta[rownames(meta),]
  meta$HTODemux = as.character(HTODemux_meta$HTODemux)
  
  MULTIseqDemux_meta<-readRDS(paste0("MULTIseqDemux/", cur_sample, ".MULTIseqDemux.rds"))@meta.data
  stopifnot(nrow(MULTIseqDemux_meta) == nrow(meta))
  MULTIseqDemux_meta<-MULTIseqDemux_meta[rownames(meta),]
  meta$MULTIseqDemux = as.character(MULTIseqDemux_meta$MULTIseqDemux)
  
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
  meta$GMM_Demux = as.character(gmm$GMM_demux)

  bff_raw_meta<-readRDS(paste0("bff_raw/", cur_sample, ".bff_raw.rds"))@meta.data
  bff_raw_meta$BFF_raw <- gsub("_", "-", bff_raw_meta$bff_raw)
  stopifnot(nrow(bff_raw_meta) == nrow(meta))
  bff_raw_meta<-bff_raw_meta[rownames(meta),]
  meta$BFF_raw = as.character(bff_raw_meta$BFF_raw)
  
  bff_cluster_file = paste0("bff_cluster/", cur_sample, ".bff_cluster.rds")
  bff_cluster_meta<-readRDS(bff_cluster_file)@meta.data
  bff_cluster_meta$BFF_cluster <- gsub("_", "-", bff_cluster_meta$bff_cluster)
  stopifnot(nrow(bff_cluster_meta) == nrow(meta))
  bff_cluster_meta<-bff_cluster_meta[rownames(meta),]
  meta$BFF_cluster = as.character(bff_cluster_meta$BFF_cluster)

  demuxmix_meta<-readRDS(paste0("demuxmix/", cur_sample, ".demuxmix.rds"))@meta.data
  stopifnot(nrow(demuxmix_meta) == nrow(meta))
  demuxmix_meta<-demuxmix_meta[rownames(meta),]
  meta$demuxmix = as.character(demuxmix_meta$demuxmix)

  hashedDrops_meta<-readRDS(paste0("hashedDrops/", cur_sample, ".hashedDrops.rds"))@meta.data
  stopifnot(nrow(hashedDrops_meta) == nrow(meta))
  hashedDrops_meta<-hashedDrops_meta[rownames(meta),]
  meta$hashedDrops = as.character(hashedDrops_meta$hashedDrops)
  
  write.csv(meta, file = paste0(cur_sample, ".results.csv"), row.names = T)
  
  obj@meta.data<-meta

  final_file = paste0(cur_sample, ".results_obj.rds")
  saveRDS(obj, final_file)
}
