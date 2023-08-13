library(Seurat)
library(ggplot2)
library(patchwork)
library(kableExtra)
library(dplyr)
library(mclust)
library(stringr)
library(aricode)

source("/home/shengq2/program/scDemultiplex_analysis/common.r")

get_plot<-function(glist1, nolegend=FALSE){
  if(nolegend){
    g1=wrap_plots(glist1, ncol=4)
  }else{
    g1=wrap_plots(glist1, ncol=4) & theme(legend.position = "right")
    g1=g1 + plot_layout(guides = "collect")
  }
  return(g1)
}

refine_cols<-names(htonames)
refine_cols<-refine_cols[!(refine_cols %in% c("scDemultiplex", "genetic_HTO"))]

prepare_sample<-function(root_dir, hashtag_to_truth, refine_cols, name){
  obj_rds_file<-paste0(root_dir, name, "/scDemultiplex_refine/", name, ".scDemultiplex_refine.rds")
  htos<-readRDS(obj_rds_file)

  meta = htos@meta.data
  truth_column = get_ground_truth_column(htos)

  cur_cols = c()
  col=refine_cols[1]
  for(col in refine_cols){
    refine_col = paste0(col, "_scDemultiplex")
    if(refine_col %in% colnames(htos@meta.data)){
      cur_cols=c(cur_cols, col, refine_col)
    }
  }

  res = check_performance(name, meta, truth_column, hashtag_to_truth, cur_cols, allow_call_name_missing = FALSE)

  return(list(htos=htos, refine_cols=cur_cols, res=res))
}

keep_digits<-function(DF, digits=3){
  is.num <- sapply(DF, is.numeric)
  DF[is.num] <- lapply(DF[is.num], round, digits)
  return(DF)
}

to_df<-function(ari_df){
  oari_df = ari_df[!grepl("scDemultiplex", rownames(ari_df)),,drop=F]
  sari_df = ari_df[grepl("scDemultiplex", rownames(ari_df)),,drop=F]
  ari_df = cbind(oari_df, sari_df)
  colnames(ari_df) = c("original", "refine")
  ari_df = data.frame(round(ari_df, digits=3))
  ari_df$combined = paste0(ari_df[,"refine"], "(", ari_df[,"original"], ")")
  return(ari_df)
}
  
name="batch1_c1"
process_sample<-function(root_dir, hashtag_to_truth, refine_cols, name){
  setwd(file.path(root_dir, name, "scDemultiplex_refine"))

  height=3000

  cat("\n\n##", name, "\n\n")
  hto_list<-prepare_sample(root_dir, hashtag_to_truth, refine_cols, name)
  obj<-hto_list$htos
  meta<-obj@meta.data
  refine_cols = hto_list$refine_cols
  res = hto_list$res

  ari_df = to_df(res$ari_df)
  write.csv(ari_df, paste0(name, ".refine.ari.csv"))

  fscore_df = to_df(res$fscore_df)
  write.csv(ari_df, paste0(name, ".refine.fscore.csv"))
  
  return(list(ari_df = ari_df, fscore_df = fscore_df))
}

all_ari = NULL
all_fscore = NULL
name=samples[1]
for(name in samples){
  lst = process_sample(root_dir, hashtag_to_truth, refine_cols, name)

  cur_ari = lst$ari_df[,"combined",drop=F]
  colnames(cur_ari) = name
  
  cur_fscore = lst$fscore_df[,"combined",drop=F]
  colnames(cur_fscore) = name
  
  if(is.null(all_ari)){
    all_ari = cur_ari
    all_fscore = cur_fscore
  }else{
    all_ari = cbind(all_ari, cur_ari[rownames(all_ari),,drop=F])
    all_fscore = cbind(all_fscore, cur_fscore[rownames(all_fscore),,drop=F])
  }
}

write.csv(all_ari, paste0(root_dir, "refine_ari.csv"), na="")
write.csv(all_fscore, paste0(root_dir, "refine_fscore.csv"), na="")
