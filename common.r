load_install<-function(library_name, library_sources=library_name){
  if(!require(library_name, character.only = T)){
    BiocManager::install(library_sources, ask=FALSE)
  }
  library(library_name, character.only = T)
}

load_install("zoo")
load_install("R.utils")
load_install("reshape2")
load_install("Matrix")
load_install("data.table")
load_install("dplyr")

load_install("Seurat")
load_install("demuxmix")
load_install("aricode")
load_install("tictoc")
load_install("cellhashR", 'BimberLab/cellhashR')
load_install("scDemultiplex", c('shengqh/cutoff', 'shengqh/scDemultiplex'))

is_unix=.Platform$OS.type == "unix"
if(is_unix) {
  root_dir="/nobackup/h_cqs/collaboration/20230522_scrna_hto/"
} else {
  root_dir="C:/projects/nobackup/h_cqs/collaboration/20230522_scrna_hto/"
}
setwd(root_dir)

scDemultiplex.p.cuts = c(0.001)
names(scDemultiplex.p.cuts) = c("scDemultiplex")

samples = c(
  "batch1_c1", 
  "batch1_c2", 
  "batch2_c1",
  "batch2_c2",
  "batch3_c1",
  "batch3_c2"
  )

sample_tags = list(
  "batch1_c1" = "Human-HTO-1,Human-HTO-2,Human-HTO-3,Human-HTO-4,Human-HTO-5,Human-HTO-6,Human-HTO-7,Human-HTO-8",
  "batch1_c2" = "Human-HTO-1,Human-HTO-2,Human-HTO-3,Human-HTO-4,Human-HTO-5,Human-HTO-6,Human-HTO-7,Human-HTO-8",
  "batch2_c1" = "Human-HTO-6,Human-HTO-7,Human-HTO-9,Human-HTO-10,Human-HTO-12,Human-HTO-13,Human-HTO-14,Human-HTO-15",
  "batch2_c2" = "Human-HTO-6,Human-HTO-7,Human-HTO-9,Human-HTO-10,Human-HTO-12,Human-HTO-13,Human-HTO-14,Human-HTO-15",
  "batch3_c1" = "Human-HTO-6,Human-HTO-7,Human-HTO-9,Human-HTO-10,Human-HTO-12,Human-HTO-13,Human-HTO-14,Human-HTO-15",
  "batch3_c2" = "Human-HTO-6,Human-HTO-7,Human-HTO-9,Human-HTO-10,Human-HTO-12,Human-HTO-13,Human-HTO-14,Human-HTO-15"
)

if(is_unix) {
  htonames<-c(
    "scDemultiplex"="scDemultiplex", 
    "HTODemux"="HTODemux", 
    "MULTIseqDemux"="MULTIseqDemux", 
    "GMM_Demux"="GMM-Demux", 
    "BFF_raw"="BFF_raw", 
    "BFF_cluster"="BFF_cluster", 
    "demuxmix"="demuxmix", 
    "hashedDrops"="hashedDrops")
  htocols=names(htonames)
  htonames[["genetic_HTO"]] = "genetic_HTO"
}else{
  htonames<-c(
    "scDemultiplex"="scDemultiplex"
  )
  htocols=names(htonames)
}

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
  htos = htos[!(htos %in% c("Negative", "Doublet", "Multiplet"))]

  fscores = unlist(lapply(htos, function(HTO){
    calculate_fscore_HTO(HTO, genetic_HTO, calls)
  }))

  return(mean(fscores))
}

get_scDemultiplex_folder<-function(root_dir, cur_sample, p.cut){
  return(paste0(root_dir, cur_sample, "/scDemultiplex.", p.cut))
}
