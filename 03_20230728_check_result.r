rm(list=ls()) 

library(Seurat)
library(ggplot2)
library(patchwork)
library(kableExtra)
library(dplyr)
library(mclust)
library(stringr)
library(scales)
library(RColorBrewer)

is_unix=.Platform$OS.type == "unix"
if(is_unix) {
  source('/home/shengq2/program/scDemultiplex_analysis/common.r')
} else {
  source('C:/Users/sheng/Programs/scDemultiplex_analysis/common.r')
}

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

get_hue_colors<-function(n, random_colors=TRUE, random.seed=20220606){
  ccolors<-hue_pal()(n)
  if(random_colors){
    x <- .Random.seed
    set.seed(random.seed)
    ccolors<-sample(ccolors, size=n)
    .Random.seed <- x
  }
  return(ccolors)
}

get_all_hto_methods<-function(root_dir, name, htocols){
  rds_file=paste0(root_dir, name, "/", name, ".results_obj.rds")
  htos<-readRDS(rds_file)
  #all_htocols = c("ground_truth", htocols)
  all_htocols = htocols
  valid_htocols = all_htocols[all_htocols %in% colnames(htos@meta.data)]
  return(valid_htocols)
}

get_alltags_umap<-function(htos, group.by, colors, legend_type, group.title=group.by){
  g = DimPlot(htos, reduction="umap", group.by=group.by) + 
    ggtitle(group.title) + 
    scale_color_manual(values=colors) + 
    theme(aspect.ratio = 1)

  if(legend_type == "corner"){
    g=g+theme(legend.position = c(2.8, -2.3))
  }else if(legend_type == "none"){
    g=g+NoLegend()
  }

  return(g)
}

get_global_name<-function(group.by){
  newcol.global.name<-gsub("-", "_", paste0(group.by, ".global"))  
  return(newcol.global.name)
}

get_global_umap<-function(htos, group.by, global_colors, legend_type, group.title=group.by){
  newcol.global.name<-get_global_name(group.by)
  newcol.global<-htos@meta.data[,group.by]
  newcol.global[!newcol.global %in% c("Negative", "Multiplet")]<-"Singlet"
  table(newcol.global)
  htos<-AddMetaData(htos, newcol.global, col.name = newcol.global.name)
  
  celllist<-list()
  celllist$Negative<-colnames(htos)[newcol.global == "Negative"]
  celllist$Multiplet<-colnames(htos)[newcol.global == "Multiplet"]
  
  g<-DimPlot(htos, reduction="umap", cells.highlight = celllist, group.by=newcol.global.name) +
    ggtitle(group.title)
  g$data$highlight=dplyr::recode(g$data$highlight, Unselected="Singlet")
  g<-g+scale_color_manual(values=global_colors)

  if(legend_type == "corner"){
    g=g+theme(legend.position = c(2.8, -2.3))
  }else if(legend_type == "none"){
    g=g+NoLegend()
  }

  return(g)
}

name="pbmc8"
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
  
  levels = c(tag_names, "Negative", "Multiplet")

  truth_column = get_ground_truth_column(htos)
  
  #replace all Doublet as Multiplex
  col=all_cols[1]
  for(col in c(all_cols, truth_column)){
    htos@meta.data[,col]<-as.character(htos@meta.data[,col])
    htos@meta.data[,col][htos@meta.data[,col]=="Doublet"]<-"Multiplet"
    #htos@meta.data[,col]<-factor(htos@meta.data[,col], levels=levels)
  }

  meta<-htos@meta.data

  colors = get_hue_colors(length(levels))
  names(colors) = levels
  global_colors=c("Singlet"="gray", "Negative"="blue", "Multiplet"="red")
  
  col = cur_htocols[1]
  for(col in cur_htocols){
    htos@meta.data[,col] = as.character(htos@meta.data[,col])

    g1 = get_alltags_umap(htos = htos, 
      group.by = col, 
      colors = colors, 
      legend_type = ifelse(col == cur_htocols[1], "corner", "none"),
      group.title = htonames[col])
    glist1[[col]] = g1

    g2 = get_global_umap(htos = htos, 
      group.by = col, 
      global_colors = global_colors, 
      legend_type = ifelse(col == cur_htocols[1], "corner", "none"),
      group.title = htonames[col])

    newcol.global.name<-get_global_name(col)
    glist2[[newcol.global.name]] = g2
  }
  
  g1<-get_plot(glist1, nolegend = TRUE)
  # png(paste0(name, ".demulti1.png"), width=4000, height=3000, res=300)
  # print(g1)
  # dev.off()

  g2<-get_plot(glist2, nolegend = TRUE)

  return(list(htos=htos, htocols=cur_htocols, g1=g1, g2=g2, colors=colors, global_colors=global_colors))
}

#name="barnyard"
#name="pbmc8"
name="batch1_c1"
process_sample<-function(root_dir, sample_tags, hashtag_to_truth, name){
  setwd(file.path(root_dir, name))

  raw_exps_rds = paste0(name, "_umi_mtx.rds")
  if(file.exists(raw_exps_rds)){
    raw_htos = readRDS(paste0(name, "_hto_mtx.rds"))
    raw_exps = readRDS(raw_exps_rds)
    common_cells = intersect(colnames(raw_htos), colnames(raw_exps))
    n_hto = ncol(raw_htos)
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
  htos<-hto_list$htos
  meta<-htos@meta.data
  meta$log_nCount_HTO<-log2(meta$nCount_HTO + 1)
  colors<-hto_list$colors
  global_colors<-hto_list$global_colors

  truth_column = get_ground_truth_column(htos)

  cat("\n\n### HTO demultiplex methods\n\n")

  cur_htocols = hto_list$htocols 

  alltb<-NULL
  col = cur_htocols[1]
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
  colnames(alltb)<-gsub("Human-", "", colnames(alltb))
  write.csv(alltb, paste0(name, ".cell.csv"))
  
  width=3000
  height=3000

  png(paste0(name, ".demulti1.png"), width=width, height=height, res=300)
  print(hto_list$g1)
  dev.off()

  png(paste0(name, ".demulti2.png"), width=width, height=height, res=300)
  print(hto_list$g2)
  dev.off()

  all_names = hto_list$htocols
  all_names = all_names[all_names != truth_column]

  #ground truth
  if(name == "barnyard"){
    truth_colors = list("HEK"="red", "MEF"="blue", "Jurkats"="purple")
  }else{
    truth_colors = colors;

    g2 = get_global_umap(htos = htos, 
      group.by = truth_column, 
      global_colors = global_colors, 
      legend_type = "default", 
      group.title =  "Ground truth")

    png(paste0(name, ".ground_truth.2.png"), width=1500, height=1300, res=300)
    print(g2)
    dev.off()
  }

  g1 = get_alltags_umap(htos = htos, 
    group.by = truth_column, 
    colors = truth_colors, 
    legend_type = "default", 
    group.title = "Ground truth")

  png(paste0(name, ".ground_truth.1.png"), width=1500, height=1300, res=300)
  print(g1)
  dev.off()

  iterations = paste0('X', 0:10)
  call_names = c(all_names, iterations)

  res = check_performance(name, meta, truth_column, hashtag_to_truth, call_names, allow_call_name_missing = TRUE)

  ari_df=res$ari_df[all_names,,drop=F]
  fscore_df=res$fscore_df[all_names,,drop=F]
  fscores_df=res$fscores_df[res$fscores_df$method %in% all_names,,drop=F]

  write.csv(ari_df, paste0(name, ".ari.csv"))
  write.csv(fscore_df, paste0(name, ".fscore.csv"))
  write.csv(fscores_df, paste0(name, ".fscore_detail.csv"), row.names=F)
  
  g<-ggplot(fscores_df, aes(method, fscore)) + geom_violin() + geom_jitter(width=0.1) + xlab("") + ylab("F score")+ ggtitle(name) + theme_bw3(rotatex=TRUE)
  png(paste0(name, ".fscore_detail.png"), width=1500, height=1300, res=300)
  print(g)
  dev.off()

  fscores_df$method = factor(fscores_df$method, levels=unique(fscores_df$method))
  fscore_long_df = reshape2::acast(fscores_df, method ~ hashtag, value.var="fscore")
  write.csv(fscore_long_df, paste0(name, ".fscore_detail_long.csv"), row.names=T)

  cur_iterations = iterations[iterations %in% rownames(res$ari_df)]
  cur_iterations = cur_iterations[1:(length(cur_iterations)-1)]
  ari_df=res$ari_df[cur_iterations,,drop=F]
  fscore_df=res$fscore_df[cur_iterations,,drop=F]
  fscores_df=res$fscores_df[res$fscores_df$method %in% cur_iterations,,drop=F]

  write.csv(ari_df, paste0(name, ".iteration.ari.csv"))
  write.csv(fscore_df, paste0(name, ".iteration.fscore.csv"))
  write.csv(fscores_df, paste0(name, ".iteration.fscore_detail.csv"), row.names=F)

  fscores_df$method = as.numeric(gsub("X", "", fscores_df$method))
  g<-ggplot(fscores_df, aes(method, fscore)) + geom_point(aes(color=hashtag)) + geom_path(aes(color=hashtag)) + xlab("") + ylab("F score")+ ggtitle(name) + theme_bw3(rotatex=TRUE)
  png(paste0(name, ".iteration.fscore_detail.png"), width=1500, height=1300, res=300)
  print(g)
  dev.off()
}

name = "pbmc8"
for(name in samples){
  process_sample(root_dir, sample_tags, hashtag_to_truth, name)
}

ari_df = NULL
fscore_df = NULL
fscore_detail_df = NULL

name = "pbmc8"
for(name in samples){
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

  fscore_detail_file = paste0(root_dir, name, "/", name, ".fscore_detail.csv")
  if(file.exists(fscore_file)){
    cur_fscore = read.csv(fscore_file, row.names=1)
    cur_fscore$dataset = name
    fscore_detail_df = rbind(fscore_detail_df, cur_fscore)
  }
}

write.csv(ari_df, paste0(root_dir, "ari.csv"), row.names=TRUE)
write.csv(fscore_df, paste0(root_dir, "fscore.csv"), row.names=TRUE)
