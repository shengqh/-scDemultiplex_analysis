library(scDemultiplex)
library(zoo)
library("R.utils")
library(reshape2)
library(Matrix)
library(data.table)
library(Seurat)
library(demuxmix)
library(tictoc)
library(dplyr)

root_dir="/nobackup/h_cqs/collaboration/20230301_scrna_hto/"
setwd(root_dir)

#samples = c("hto12", "pbmc", "batch1", "batch2", "batch3")
samples = c(
  "hto12", 
  "pbmc", 
  "batch1", 
  "batch2", 
  "batch3")

sample_tags = list(
  "hto12" = "HEK_A,HEK_B,HEK_C,K562_A,K562_B,K562_C,KG1_A,KG1_B,KG1_C,THP1_A,THP1_B,THP1_C",
  "pbmc" = "HTO_A,HTO_B,HTO_C,HTO_D,HTO_E,HTO_F,HTO_G,HTO_H",
  "batch1" = "Human-HTO-1,Human-HTO-2,Human-HTO-3,Human-HTO-4,Human-HTO-5,Human-HTO-6,Human-HTO-7,Human-HTO-8",
  "batch2" = "Human-HTO-6,Human-HTO-7,Human-HTO-9,Human-HTO-10,Human-HTO-12,Human-HTO-13,Human-HTO-14,Human-HTO-15",
  "batch3" = "Human-HTO-6,Human-HTO-7,Human-HTO-9,Human-HTO-10,Human-HTO-12,Human-HTO-13,Human-HTO-14,Human-HTO-15"
)

htonames<-c(
  "scDemultiplex_cutoff"="scDemultiplex_cutoff", 
  "scDemultiplex"="scDemultiplex", 
  "HTODemux"="HTODemux", 
  "MULTIseqDemux"="MULTIseqDemux", 
  "GMM_Demux"="GMM-Demux", 
  "BFF_raw"="BFF_raw", 
  #"BFF_cluster"="BFF_cluster", 
  "demuxmix"="demuxmix", 
  "hashedDrops"="hashedDrops")
htocols=names(htonames)
htonames[["genetic_HTO"]] = "genetic_HTO"

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

calculate_fscore_HTO<-function(HTO, genetic_HTO, calls){
  tp <- sum(calls == HTO & genetic_HTO == HTO) #True positive rate
  fp <- sum(calls == HTO & genetic_HTO != HTO) #False positive rate
  fn <- sum(calls != HTO & genetic_HTO == HTO) #False negative rate
  f <- tp / (tp + 0.5 * (fp + fn))
  return(f)
}

calculate_fscore<-function(genetic_HTO, calls){
  genetic_HTO = as.character(genetic_HTO)
  calls = as.character(calls)

  htos = unique(genetic_HTO)
  htos = htos[!(htos %in% c("Negative", "Doublet", "Multiplex"))]
  
  fscores = unlist(lapply(htos, function(HTO){
    calculate_fscore_HTO(HTO, genetic_HTO, calls)
  }))

  return(mean(fscores))
}
