library(Seurat)
library(ggplot2)
library(patchwork)
library(kableExtra)
library(dplyr)
library(mclust)
library(stringr)
library(aricode)

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
    g1=wrap_plots(glist1, ncol=4)
  }else{
    g1=wrap_plots(glist1, ncol=4) & theme(legend.position = "right")
    g1=g1 + plot_layout(guides = "collect")
  }
  return(g1)
}

refine_cols<-names(htonames)
refine_cols<-refine_cols[!(refine_cols %in% c("scDemultiplex", "genetic_HTO"))]

prepare_sample<-function(root_dir, sample_tags, name){
  obj_rds_file<-paste0(root_dir, name, "/scDemultiplex_refine/", name, ".scDemultiplex_refine.rds")
  htos<-readRDS(obj_rds_file)

  glist1<-list()
  glist2<-list()

  cur_cols = c("genetic_HTO", "scDemultiplex")
  col=refine_cols[1]
  for(col in refine_cols){
    refine_col = paste0(col, "_scDemultiplex")
    if(refine_col %in% colnames(htos@meta.data)){
      cur_cols=c(cur_cols, col, refine_col)
    }
  }

  tag_names=gsub('_', '-', unlist(str_split(sample_tags[[name]], ',')))

  #replace all Doublet as Multiplex
  for(col in cur_cols){
    htos@meta.data[,col]<-as.character(htos@meta.data[,col])
    htos@meta.data[,col][htos@meta.data[,col]=="Doublet"]<-"Multiplet"
    htos@meta.data[,col]<-factor(htos@meta.data[,col], levels=c(tag_names, "Negative", "Multiplet"))
  }

  meta<-htos@meta.data

  for(col in cur_cols){
    newcol.global.name<-gsub("-", "_", paste0(col, ".global"))
    newcol.global<-as.character(meta[,col])
    newcol.global[!newcol.global %in% c("Negative", "Multiplet")]<-"Singlet"
    htos<-AddMetaData(htos, newcol.global, col.name = newcol.global.name)

    glist1[[col]] = DimPlot(htos, reduction="umap", group.by=col) + ggtitle(col)
    
    celllist<-list()
    celllist$Negative<-colnames(htos)[newcol.global == "Negative"]
    celllist$Multiplex<-colnames(htos)[newcol.global == "Multiplet"]
    
    glist2[[newcol.global.name]] = DimPlot(htos, reduction="umap", cells.highlight = celllist, group.by=newcol.global.name) + 
      ggtitle(col) + 
      scale_color_manual(labels=c('Singlet', 'Negative', 'Multiplet'), values =c("gray", "blue", "red"))
  }
  
  g1<-get_plot(glist1)
  g2<-get_plot(glist2)

  return(list(htos=htos, refine_cols=cur_cols, g1=g1, g2=g2))
}

name="batch1_c1"
process_sample<-function(root_dir, sample_tags, name){
  setwd(file.path(root_dir, name, "scDemultiplex_refine"))

  height=3000

  cat("\n\n##", name, "\n\n")
  hto_list<-prepare_sample(root_dir, sample_tags, name)
  obj<-hto_list$htos
  meta<-obj@meta.data
  meta$log_nCount_HTO<-log2(meta$nCount_HTO + 1)

  refine_cols = hto_list$refine_cols

  cat("\n\n### HTO demultiplex methods\n\n")
  
  alltb<-NULL
  col=hto_list$refine_cols[1]
  for (col in hto_list$refine_cols){
    coltb<-table(meta[,col])
    if(is.null(alltb)){
      alltb<-rbind(alltb, coltb)
    }else{
      coltb<-coltb[colnames(alltb)]
      alltb<-rbind(alltb, coltb)
    }
  }
  rownames(alltb)<-hto_list$refine_cols
  alltb<-t(alltb)
  write.csv(alltb, paste0(name, ".cell.csv"))
  
  png(paste0(name, ".demulti1.png"), width=4300, height=height, res=300)
  print(hto_list$g1)
  dev.off()

  png(paste0(name, ".demulti2.png"), width=4300, height=height, res=300)
  print(hto_list$g2)
  dev.off()

  methods = refine_cols[refine_cols != "genetic_HTO"]
  aris = unlist(lapply(methods, function(x){
    ARI(meta[,"genetic_HTO"], meta[,x])
  }))
  fscores = unlist(lapply(methods, function(x){
    calculate_fscore(meta[,"genetic_HTO"], meta[,x])
  }))
  df=data.frame(ARI=aris, Fscore=fscores)
  rownames(df) = methods
  write.csv(df, paste0(name, ".score.csv"))
  
  return(df)
}

all_ari = NULL
all_fscore = NULL
name=samples[1]
for(name in samples){
  df = process_sample(root_dir, sample_tags, name)

  cur_ari = df[,"ARI",drop=F]
  colnames(cur_ari) = name
  
  cur_fscore = df[,"Fscore",drop=F]
  colnames(cur_fscore) = name
  
  if(is.null(all_ari)){
    all_ari = cur_ari
    all_fscore = cur_fscore
  }else{
    all_ari = cbind(all_ari, cur_ari[rownames(all_ari),,drop=F])
    all_fscore = cbind(all_fscore, cur_fscore[rownames(all_fscore),,drop=F])
  }
}

write.csv(all_ari, paste0(root_dir, "refine_ari.csv"))
write.csv(all_fscore, paste0(root_dir, "refine_fscore.csv"))
