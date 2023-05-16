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

prepare_sample<-function(root_dir, sample_tags, name){
  rds_file=paste0(root_dir, name, "/", name, ".results_obj.rds")

  htos<-readRDS(rds_file)

  glist1<-list()
  glist2<-list()

  col=htocols[1]

  tag_names=gsub('_', '-', unlist(str_split(sample_tags[[name]], ',')))

  if(grepl("batch", name)){
    cur_htocols = c("genetic_HTO", htocols)
  }else{
    cur_htocols = htocols
  }

  #replace all Doublet as Multiplex
  for(col in cur_htocols){
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
    
    glist2[[newcol.global.name]] = DimPlot(htos, reduction="umap", cells.highlight = celllist, group.by=newcol.global.name) + 
      ggtitle(htonames[col]) + 
      scale_color_manual(labels=c('Singlet', 'Negative', 'Multiplet'), values =c("gray", "blue", "red"))
  }
  
  g1<-get_plot(glist1)
  g2<-get_plot(glist2)

  return(list(htos=htos, htocols=htocols, g1=g1, g2=g2))
}

name="batch1"
process_sample<-function(root_dir, sample_tags, name){
  setwd(file.path(root_dir, name))

  if(!grepl("batch", name)){
    raw_htos = readRDS(paste0(name, "_hto_mtx.rds"))
    raw_exps = readRDS(paste0(name, "_umi_mtx.rds"))
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

  if(grepl("batch", name)){
    cur_htocols = c("genetic_HTO", htocols)
  }else{
    cur_htocols = htocols
  }

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
  
  cat("\n\n### DE analysis\n\n")
  alltb<-NULL
  col<-hto_list$htocols[1]
  for (col in hto_list$htocols){
    print(col)
    subobj<-obj
    Idents(subobj)<-col
    markers<-FindAllMarkers(subobj, pseudocount.use = FALSE,only.pos = TRUE)
    markers<-markers[!(markers$cluster %in% c("Negative", "Multiplet")),]
    markers<-markers[order(markers$avg_log2FC, decreasing=TRUE),c("cluster", "gene", "avg_log2FC"), drop=F]
    markers<-markers[!duplicated(markers$cluster),,drop=F]
    markers<-markers[order(markers$cluster),]
    rownames(markers)<-markers$cluster
    markers<-markers[,"avg_log2FC", drop=F]
    if(is.null(alltb)){
      alltb<-markers
    }else{
      alltb<-cbind(alltb, markers[rownames(alltb),])
    }
  }
  colnames(alltb)<-hto_list$htocols
  write.csv(alltb, paste0(name, ".de.csv"))
  
  if(grepl("batch", name)){
    ari_list=list()
    col=hto_list$htocols[1]
    for(col in hto_list$htocols){
      ari=adjustedRandIndex(meta[,col], meta[,"genetic_HTO"])
      ari_list[[col]] = ari
    }
    ari_df<-t(data.frame(ari_list))
    colnames(ari_df)<-"ari"
    write.csv(ari_df, paste0(name, ".ari.csv"))

    fscore_list=list()
    col=hto_list$htocols[1]
    for(col in hto_list$htocols){
      fscore=calculate_fscore(meta[,"genetic_HTO"], meta[,col])
      fscore_list[[col]] = fscore
    }
    fscore_df<-t(data.frame(fscore_list))
    colnames(fscore_df)<-"Fscore"
    write.csv(fscore_df, paste0(name, ".fscore.csv"))
  }

  if(name == "hto12"){
    hto12_combined_obj_rds<-paste0(root_dir, "/hto12/hto12.combined.rds")
    expobj<-readRDS(hto12_combined_obj_rds)
    expobj<-expobj[,colnames(obj)]
    ari_list=list()
    col=hto_list$htocols[1]
    for(col in hto_list$htocols){
      df = data.frame("exp_cluster"=expobj$seurat_clusters, "hto_assign"=unlist(FetchData(obj, col)))
      write.csv(df, paste0(name, ".", col, ".ari.csv"))
      
      ari=adjustedRandIndex(expobj$seurat_clusters, unlist(FetchData(obj, col)))
      ari_list[[col]] = ari
    }
    ari_df<-t(data.frame(ari_list))
    colnames(ari_df)<-"ari"
    write.csv(ari_df, paste0(name, ".ari.csv"))

    stopifnot(all(colnames(obj) == colnames(expobj)))

    glist=list()
    newcols=c()
    for(col in hto_list$htocols){
      newcol = paste0(col, "_celltype")
      newcolvalue = gsub("-.", "", unlist(FetchData(obj, col)))
      newcolvalue = factor(newcolvalue, levels=c("HEK", "K562", "KG1", "THP1", "Negative", "Multiplet"))

      obj<-AddMetaData(obj, newcolvalue, col.name = newcol)
      expobj<-AddMetaData(expobj, newcolvalue, col.name = newcol)

      cts<-sort(unique(newcolvalue))      
      newcols<-c(newcols, newcol)
      glist[[col]] = DimPlot(expobj, reduction = "umap", group.by=newcol) + ggtitle(htonames[col])
    }
    g1=get_plot(glist)
    png("hto12.exp_validation.all.png", width=width, height=height, res=300)
    print(g1)
    dev.off()

    ct=cts[1]
    for(ct in cts){
      glist=list()
      for(col in hto_list$htocols){
        newcol = paste0(col, "_celltype")
        newcolvalue = unlist(FetchData(obj, newcol))
        cells<-colnames(obj)[newcolvalue == ct]
        gname = paste0(col, ":", ct)
        g<-DimPlot(expobj, reduction = "umap", cells.highlight = cells, cols.highlight = "red") + ggtitle(htonames[col]) + NoLegend()
        glist[[gname]] = g
      }
      g1=get_plot(glist, nolegend=TRUE) + plot_annotation(title=ct, theme = theme(plot.title = element_text(size = 20, hjust = 0.5)))

      png(paste0("hto12.exp_validation.", ct, ".png"), width=3000, height=height, res=300)
      print(g1)
      dev.off()
    }
  }
}

for(name in samples){
  process_sample(root_dir, sample_tags, name)
}
