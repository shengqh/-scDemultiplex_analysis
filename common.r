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
load_install("aricode")
load_install("tictoc")
load_install("mclust")

is_unix=.Platform$OS.type == "unix"
if(is_unix) {
  root_dir="/nobackup/h_cqs/collaboration/20230725_scrna_hto/"
} else {
  root_dir="C:/projects/nobackup/h_cqs/collaboration/20230725_scrna_hto/"
}
if(!dir.exists(root_dir)){
  dir.create(root_dir)
}
setwd(root_dir)

scDemultiplex.p.cuts = c(0.001)
names(scDemultiplex.p.cuts) = c("scDemultiplex")

samples = c(
  "barnyard",
  "pbmc8",
  "batch1_c1", 
  "batch1_c2", 
  "batch2_c1",
  "batch2_c2",
  "batch3_c1",
  "batch3_c2"
)

batch1_tag = "BAL-01,BAL-02,BAL-03,BAL-04,BAL-05,BAL-06,BAL-07,BAL-08"
batch2_tag = "BAL-09,BAL-10,BAL-11,BAL-12,BAL-13,BAL-14,BAL-15,BAL-16"
batch3_tag = "BAL-17,BAL-18,BAL-19,BAL-20,BAL-21,BAL-22,BAL-23,BAL-24"
sample_tags = list(
  "batch1_c1" = batch1_tag,
  "batch1_c2" = batch1_tag,
  "batch2_c1" = batch2_tag,
  "batch2_c2" = batch2_tag,
  "batch3_c1" = batch3_tag,
  "batch3_c2" = batch3_tag,
  "barnyard" = "Bar1,Bar2,Bar3,Bar4,Bar5,Bar6,Bar7,Bar8,Bar9,Bar10,Bar11,Bar12",
  "pbmc8" = "HTO-A,HTO-B,HTO-C,HTO-D,HTO-E,HTO-F,HTO-G,HTO-H"
  #"pbmc8" = "BatchA,BatchB,BatchC,BatchD,BatchE,BatchF,BatchG,BatchH"
)

hashtag_to_truth = list(
  "barnyard" = list(
    "Bar1" = "HEK",
    "Bar2" = "HEK",
    "Bar3" = "MEF",
    "Bar4" = "MEF",
    "Bar5" = "Jurkats",
    "Bar6" = "Jurkats",
    "Bar7" = "Jurkats",
    "Bar8" = "Jurkats",
    "Bar9" = "Jurkats",
    "Bar10" = "Jurkats",
    "Bar11" = "Jurkats",
    "Bar12" = "Jurkats",
    "Doublet" = "Multiplet",
    "Multiplet" = "Multiplet",
    "Negative" = "Negative"
  )
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
  htonames[["ground_truth"]] = "ground_truth"
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

calculate_fscore_HTO<-function(HTO, ground_truth, calls){
  tp <- sum(calls == HTO & ground_truth == HTO) #True positive rate
  fp <- sum(calls == HTO & ground_truth != HTO) #False positive rate
  fn <- sum(calls != HTO & ground_truth == HTO) #False negative rate
  f <- tp / (tp + 0.5 * (fp + fn))
  return(f)
}

calculate_fscore<-function(ground_truth, calls){
  ground_truth = as.character(ground_truth)
  calls = as.character(calls)

  htos = unique(ground_truth)
  htos = htos[!(htos %in% c("Negative", "Doublet", "Multiplet", "unassigned", "doublet"))]

  fscores = unlist(lapply(htos, function(HTO){
    calculate_fscore_HTO(HTO, ground_truth, calls)
  }))
  names(fscores) = htos

  return(fscores)
}

get_scDemultiplex_folder<-function(root_dir, cur_sample, p.cut){
  return(paste0(root_dir, cur_sample, "/scDemultiplex.", p.cut))
}

theme_bw3 <- function(rotatex=FALSE) { 
	result = theme_bw() +
    theme(
      strip.background = element_rect(fill = NA, colour = 'black'),
      panel.border = element_rect(fill = NA, color = "black"),			
      axis.line = element_line(colour = "black", linewidth = 0.5)
    )

  if(rotatex){
    result = result + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  }
  result
}

to_ground_truth<-function(hashtag_to_truth, name, meta, col){
  if(name %in% names(hashtag_to_truth)){
    actual_map = hashtag_to_truth[[name]] 
    stopifnot(all(meta[,col] %in% names(actual_map)))
    actual_values = unlist(actual_map[as.character(meta[,col])])
  }else{
    actual_values = as.character(meta[,col])
  }
  return(actual_values)
}

get_ground_truth_column <- function(htos){
  truth_column = ifelse("ground_truth" %in% colnames(htos@meta.data), "ground_truth", "genetic_HTO")  
  stopifnot(truth_column %in% colnames(htos@meta.data))
  return(truth_column)
}

check_performance<-function(name, meta, truth_column, hashtag_to_truth, call_names, allow_call_name_missing = TRUE){
  ari_list=list()
  fscore_list=list()
  fscores_list=list()
  
  col=call_names[1]
  for(col in call_names){
    if(!(col %in% colnames(meta))){
      if(allow_call_name_missing){
        next 
      }else{
        stop("No column ", col, " in ", name)
      }
    }
    calls = to_ground_truth(hashtag_to_truth, name, meta, col)
    ground_truth = meta[,truth_column]

    common_hto = intersect(unique(calls), unique(ground_truth))
    if(length(common_hto) == 0){
      stop("No common HTO between calls and ground truth for ", col, " in ", name, "\ncalls: ", paste0(unique(calls), collapse=","), "\ntruth: ", paste0(unique(ground_truth), collapse=","))
    }

    ari=adjustedRandIndex(ground_truth, calls)
    ari_list[[col]] = ari

    fscores=calculate_fscore(ground_truth, calls)
    fscore_list[[col]] = mean(fscores)
    fscores_list[[col]] = fscores
  }

  ari_df<-t(data.frame(ari_list))
  colnames(ari_df)<-"ari"

  fscore_df<-t(data.frame(fscore_list))
  colnames(fscore_df)<-"Fscore"

  fscores_df<-NULL
  for(col in names(fscores_list)){
    fscores = fscores_list[[col]]
    cur_df = data.frame(method=col, hashtag=names(fscores), fscore=fscores)
    fscores_df = rbind(fscores_df, cur_df)
  }

  return(list(ari_df=ari_df, fscore_df=fscore_df, fscores_df=fscores_df))
}

get_date_str<-function(){
  date_str = format(Sys.time(), "%Y%m%d")
  return(date_str)
}

do_scDemultiplex<-function(root_dir, cur_sample, p.cuts=0.001, do_rlog=FALSE){
  load_install("scDemultiplex", c('shengqh/cutoff', 'shengqh/scDemultiplex'))

  #in order to run dff_cluster, we will need to install preprocessCore manually as single thread mode
  #however, it would not work for scDemultiplex. So we need to install preprocessCore as multi-thread mode again
  #devtools::install_github('bmbolstad/preprocessCore', force=TRUE)

  for(p.cut in p.cuts){
    suffix = ifelse(do_rlog, "rlog", p.cut)
    result_folder = get_scDemultiplex_folder(root_dir, cur_sample, suffix)
    if(!dir.exists(result_folder)){
      dir.create(result_folder)
    }
    setwd(result_folder)
  
    final_rds=paste0(cur_sample, ".scDemultiplex.rds")
    if(!file.exists(final_rds)){
      message("scDemultiplex", p.cut, "...\n")
  
      rds_file=paste0("../", cur_sample, ".obj.rds")
      obj=readRDS(rds_file)
      if(do_rlog){
        obj <- NormalizeData(obj, normalization.method = "LogNormalize")
      }
  
      ntags = nrow(obj)
      
      output_prefix<-paste0(cur_sample, ".HTO")
      
      message(paste0("starting ", cur_sample, " cutoff ...\n"))
      tic()
      cat("  scDemultiplex_cutoff ...\n")
      obj<-demulti_cutoff(obj, output_prefix=output_prefix, cutoff_startval = 0, mc.cores=ntags)
      toc1=toc()

      obj<-hto_plot(obj, paste0(output_prefix, ".cutoff"), group.by="scDemultiplex_cutoff")

      if(!do_rlog){
        message(paste0("starting ", cur_sample, " cutoff ...\n"))
        tic()
        obj<-demulti_refine(obj, output_prefix=output_prefix, p.cut=p.cut, refine_negative_doublet_only=FALSE, mc.cores=ntags)
        obj$scDemultiplex_full=obj$scDemultiplex
        obj$scDemultiplex_full.global=obj$scDemultiplex.global
        toc3=toc()
  
        saveRDS(list("cutoff"=toc1, "full"=toc3), paste0(cur_sample, ".scDemultiplex.tictoc.rds"))
        obj<-hto_plot(obj, paste0(output_prefix, ".full_p"), group.by="scDemultiplex_full")
      }
  
      saveRDS(obj, final_rds)
    }
  }
}
