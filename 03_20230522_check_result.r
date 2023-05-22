library(Seurat)
library(ggplot2)
library(patchwork)
library(kableExtra)
library(dplyr)
library(mclust)
library(stringr)

theme_bw3 <- function() { 
	theme_bw() +
	theme(
		strip.background = element_rect(fill = NA, colour = 'black'),
		panel.border = element_rect(fill = NA, color = "black"),			
		axis.line = element_line(colour = "black", size = 0.5)
	)
}

source("/home/shengq2/program/scDemultiplex_analysis/common.r")

get_plot<-function(glist1, nolegend=FALSE){
  if(nolegend){
    g1=wrap_plots(glist1, ncol=3)
  }else{
    g1=wrap_plots(glist1, ncol=3) & theme(legend.position = "right")
    g1=g1 + plot_layout(guides = "collect")
  }
  return(g1)
}

get_iteration_colnames<-function(meta){
  if("X1" %in% colnames(meta)){
    possible_names = paste0("X", c(0:10))
    valid_names = possible_names[possible_names %in% colnames(meta)]
    return(valid_names)
  }else{
    return(c())
  }
}

get_all_hto_methods<-function(root_dir, name, htocols){
  rds_file=paste0(root_dir, name, "/", name, ".results_obj.rds")
  htos<-readRDS(rds_file)
  all_htocols = c("genetic_HTO", htocols)
  valid_htocols = all_htocols[all_htocols %in% colnames(htos@meta.data)]
  return(valid_htocols)
}

prepare_sample<-function(root_dir, sample_tags, name){
  rds_file=paste0(root_dir, name, "/", name, ".results_obj.rds")

  htos<-readRDS(rds_file)

  glist1<-list()
  glist2<-list()

  cur_htocols = get_all_hto_methods(root_dir, name, htocols)
  
  tag_names=gsub('_', '-', unlist(str_split(sample_tags[[name]], ',')))

  if(0){
    iteration_names = get_iteration_colnames(htos@meta.data)
    all_cols = c(cur_htocols, iteration_names)
  }else{
    all_cols = cur_htocols
  }
  
  #replace all Doublet as Multiplex
  col=all_cols[1]
  for(col in all_cols){
    htos@meta.data[,col]<-as.character(htos@meta.data[,col])
    htos@meta.data[,col][htos@meta.data[,col]=="Doublet"]<-"Multiplet"
    htos@meta.data[,col]<-factor(htos@meta.data[,col], levels=c(tag_names, "Negative", "Multiplet"))
  }

  meta<-htos@meta.data

  for(col in cur_htocols){
    newcol.global.name<-gsub("-", "_", paste0(col, ".global"))
    newcol.global<-as.character(meta[,col])
    newcol.global[!newcol.global %in% c("Negative", "Multiplet")]<-"Singlet"
    htos<-AddMetaData(htos, newcol.global, col.name = newcol.global.name)

    glist1[[col]] = DimPlot(htos, reduction="umap", group.by=col) + ggtitle(htonames[col])
    
    celllist<-list()
    celllist$Negative<-colnames(htos)[newcol.global == "Negative"]
    celllist$Multiplex<-colnames(htos)[newcol.global == "Multiplet"]
    
    g<-DimPlot(htos, reduction="umap", cells.highlight = celllist, group.by=newcol.global.name) + 
      ggtitle(htonames[col])
      
    if("Negative" %in% meta[,col]){
      g<-g+scale_color_manual(labels=c('Singlet', 'Negative', 'Multiplet'), values =c("gray", "blue", "red"))
    }else{
      g<-g+scale_color_manual(labels=c('Singlet', 'Multiplet'), values =c("gray", "red"))
    }
    glist2[[newcol.global.name]] = g
  }
  
  g1<-get_plot(glist1)
  g2<-get_plot(glist2)

  return(list(htos=htos, htocols=cur_htocols, g1=g1, g2=g2))
}

name="batch1_c1"
process_sample<-function(root_dir, sample_tags, name){
  setwd(file.path(root_dir, name))

  raw_exps_rds = paste0(name, "_umi_mtx.rds")
  if(file.exists(raw_exps_rds)){
    raw_htos = readRDS(paste0(name, "_hto_mtx.rds"))
    raw_exps = readRDS(raw_exps_rds)
    if(name == "hto12"){
      common_cells = intersect(rownames(raw_htos), colnames(raw_exps))
      n_hto = nrow(raw_htos)
    }else{
      common_cells = intersect(colnames(raw_htos), colnames(raw_exps))
      n_hto = ncol(raw_htos)
    }
    n_umi = ncol(raw_exps)

    final_cells = colnames(readRDS(paste0(name, ".obj.rds")))
    write.csv(data.frame("Category" = c("HTO cells", "UMI cells", "Common cells", "Filtered cells"),
                        "Cell" = c(n_hto, n_umi, length(common_cells), length(final_cells))), paste0(name, ".ncell.csv"), row.names=F)
    rm(raw_htos)
    rm(raw_exps)
  }else{
    final_cells = colnames(readRDS(paste0(name, ".obj.rds")))
    write.csv(data.frame("Category" = c("HTO cells"),
                        "Cell" = c(length(final_cells))), paste0(name, ".ncell.csv"), row.names=F)
  }

  cat("\n\n##", name, "\n\n")
  hto_list<-prepare_sample(root_dir, sample_tags, name)
  obj<-hto_list$htos
  meta<-obj@meta.data
  meta$log_nCount_HTO<-log2(meta$nCount_HTO + 1)

  cat("\n\n### HTO demultiplex methods\n\n")

  cur_htocols = hto_list$htocols 

  alltb<-NULL
  for (col in cur_htocols){
    coltb<-table(meta[,col])
    if(is.null(alltb)){
      alltb<-rbind(alltb, coltb)
    }else{
      coltb<-coltb[colnames(alltb)]
      alltb<-rbind(alltb, coltb)
    }
  }
  rownames(alltb)<-cur_htocols
  alltb<-t(alltb)
  write.csv(alltb, paste0(name, ".cell.csv"))
  
  width=3300
  height=3000

  png(paste0(name, ".demulti1.png"), width=width, height=height, res=300)
  print(hto_list$g1)
  dev.off()

  png(paste0(name, ".demulti2.png"), width=width, height=height, res=300)
  print(hto_list$g2)
  dev.off()

  all_names = hto_list$htocols
  all_names = all_names[all_names != "genetic_HTO"]

  ari_list=list()
  col=all_names[1]
  for(col in all_names){
    ari=ARI(meta[,col], meta[,"genetic_HTO"])
    #ari=adjustedRandIndex(meta[,col], meta[,"genetic_HTO"])#equal result as ARI function
    ari_list[[col]] = ari
  }
  ari_df<-t(data.frame(ari_list))
  colnames(ari_df)<-"ari"
  write.csv(ari_df, paste0(name, ".ari.csv"))
  
  fscore_list=list()
  col=all_names[1]
  for(col in all_names){
    fscore=calculate_fscore(meta[,"genetic_HTO"], meta[,col])
    fscore_list[[col]] = fscore
  }
  fscore_df<-t(data.frame(fscore_list))
  colnames(fscore_df)<-"Fscore"
  write.csv(fscore_df, paste0(name, ".fscore.csv"))
}

ari_df = NULL
fscore_df = NULL
name = samples[1]
for(name in samples){
  process_sample(root_dir, sample_tags, name)
  
  ari_file = paste0(root_dir, name, "/", name, ".ari.csv")
  if(file.exists(ari_file)){
    cur_ari = read.csv(ari_file, row.names=1)
    colnames(cur_ari) = name
    if(is.null(ari_df)){
      ari_df = cur_ari
    }else{
      stopifnot(rownames(ari_df) == rownames(cur_ari))
      ari_df = cbind(ari_df, cur_ari)
    }
  }
  
  fscore_file = paste0(root_dir, name, "/", name, ".fscore.csv")
  if(file.exists(fscore_file)){
    cur_fscore = read.csv(fscore_file, row.names=1)
    colnames(cur_fscore) = name
    if(is.null(fscore_df)){
      fscore_df = cur_fscore
    }else{
      stopifnot(rownames(fscore_df) == rownames(cur_fscore))
      fscore_df = cbind(fscore_df, cur_fscore)
    }
  }
}

write.csv(ari_df, paste0(root_dir, "ari.csv"), row.names=TRUE)
write.csv(fscore_df, paste0(root_dir, "fscore.csv"), row.names=TRUE)
