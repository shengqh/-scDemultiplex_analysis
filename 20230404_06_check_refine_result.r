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

refine_cols<-c("scDemultiplex_full", paste0(c("HTODemux", "MULTIseqDemux", "GMM_Demux", "BFF_raw", "BFF_cluster"), "_scDemultiplex"))
refine_names<-c("scDemultiplex", "HTODemux", "MULTIseqDemux", "GMM-Demux", "BFF_raw", "BFF_cluster")
names(refine_names)<-refine_cols

prepare_sample<-function(root_dir, sample_tags, name){
  obj_rds_file<-paste0(root_dir, name, "/", name, ".results_obj.rds")
  htos<-readRDS(obj_rds_file)

  meta_rds_file=paste0(root_dir, name, "/scDemultiplex_refine/", name, ".BFF_cluster_scDemultiplex_p0.001.rds")
  htos@meta.data<-readRDS(meta_rds_file)

  glist1<-list()
  glist2<-list()

  col=refine_cols[1]

  tag_names=gsub('_', '-', unlist(str_split(sample_tags[[name]], ',')))

  #replace all Doublet as Multiplex
  for(col in refine_cols){
    htos@meta.data[,col]<-as.character(htos@meta.data[,col])
    htos@meta.data[,col][htos@meta.data[,col]=="Doublet"]<-"Multiplet"
    htos@meta.data[,col]<-factor(htos@meta.data[,col], levels=c(tag_names, "Negative", "Multiplet"))
  }

  meta<-htos@meta.data

  for(col in refine_cols){
    newcol.global.name<-gsub("-", "_", paste0(col, ".global"))
    newcol.global<-as.character(meta[,col])
    newcol.global[!newcol.global %in% c("Negative", "Multiplet")]<-"Singlet"
    htos<-AddMetaData(htos, newcol.global, col.name = newcol.global.name)

    glist1[[col]] = DimPlot(htos, reduction="umap", group.by=col) + ggtitle(refine_names[col])
    
    celllist<-list()
    celllist$Negative<-colnames(htos)[newcol.global == "Negative"]
    celllist$Multiplex<-colnames(htos)[newcol.global == "Multiplet"]
    
    glist2[[newcol.global.name]] = DimPlot(htos, reduction="umap", cells.highlight = celllist, group.by=newcol.global.name) + 
      ggtitle(refine_names[col]) + 
      scale_color_manual(labels=c('Singlet', 'Negative', 'Multiplet'), values =c("gray", "blue", "red"))
  }
  
  g1<-get_plot(glist1)
  g2<-get_plot(glist2)

  return(list(htos=htos, refine_cols=refine_cols, g1=g1, g2=g2))
}

name="hto12"
process_sample<-function(root_dir, sample_tags, name){
  setwd(file.path(root_dir, name, "scDemultiplex_refine"))

  cat("\n\n##", name, "\n\n")
  hto_list<-prepare_sample(root_dir, sample_tags, name)
  obj<-hto_list$htos
  meta<-obj@meta.data
  meta$log_nCount_HTO<-log2(meta$nCount_HTO + 1)

  cat("\n\n### HTO demultiplex methods\n\n")
  
  alltb<-NULL
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
  
  png(paste0(name, ".demulti1.png"), width=3300, height=2000, res=300)
  print(hto_list$g1)
  dev.off()

  png(paste0(name, ".demulti2.png"), width=3300, height=2000, res=300)
  print(hto_list$g2)
  dev.off()
  
  cat("\n\n### DE analysis\n\n")
  alltb<-NULL
  col<-hto_list$refine_cols[1]
  for (col in hto_list$refine_cols){
    print(col)
    subobj<-obj
    Idents(subobj)<-col
    markers<-FindAllMarkers(subobj, pseudocount.use = FALSE,only.pos = TRUE)
    markers<-markers[!(markers$cluster %in% c("Negative", "Multiplet")),]
    markers<-markers[order(markers$gene),c("gene", "avg_log2FC"), drop=F]
    rownames(markers)<-markers$gene
    markers<-markers[,"avg_log2FC", drop=F]
    colnames(markers)=col
    if(is.null(alltb)){
      alltb<-markers
    }else{
      alltb<-cbind(alltb, markers)
    }
  }
  write.csv(alltb, paste0(name, ".de.csv"))
  
  if(name == "hto12"){
    hto12_combined_obj_rds<-paste0(root_dir, "/hto12/hto12.combined.rds")
    expobj<-readRDS(hto12_combined_obj_rds)
    expobj<-expobj[,colnames(obj)]
    ari_list=list()
    col=hto_list$refine_cols[1]
    for(col in hto_list$refine_cols){
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
    for(col in hto_list$refine_cols){
      newcol = paste0(col, "_celltype")
      newcolvalue = gsub("-.", "", unlist(FetchData(obj, col)))
      newcolvalue = factor(newcolvalue, levels=c("HEK", "K562", "KG1", "THP1", "Negative", "Multiplet"))

      obj<-AddMetaData(obj, newcolvalue, col.name = newcol)
      expobj<-AddMetaData(expobj, newcolvalue, col.name = newcol)

      cts<-sort(unique(newcolvalue))      
      newcols<-c(newcols, newcol)
      glist[[col]] = DimPlot(expobj, reduction = "umap", group.by=newcol) + ggtitle(refine_names[col])
    }
    g1=get_plot(glist)
    png("hto12.exp_validation.all.png", width=3300, height=2000, res=300)
    print(g1)
    dev.off()

    ct=cts[1]
    for(ct in cts){
      glist=list()
      for(col in hto_list$refine_cols){
        newcol = paste0(col, "_celltype")
        newcolvalue = unlist(FetchData(obj, newcol))
        cells<-colnames(obj)[newcolvalue == ct]
        gname = paste0(col, ":", ct)
        g<-DimPlot(expobj, reduction = "umap", cells.highlight = cells, cols.highlight = "red") + ggtitle(refine_names[col]) + NoLegend()
        glist[[gname]] = g
      }
      g1=get_plot(glist, nolegend=TRUE) + plot_annotation(title=ct, theme = theme(plot.title = element_text(size = 20, hjust = 0.5)))

      png(paste0("hto12.exp_validation.", ct, ".png"), width=3000, height=2000, res=300)
      print(g1)
      dev.off()
    }
  }
}

for(name in samples){
  process_sample(root_dir, sample_tags, name)
}
