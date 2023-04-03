source('/home/shengq2/program/scDemultiplex_analysis/common.r')

hto12_combined_obj_rds<-paste0(root_dir, "/hto12/hto12.combined.rds")
hto.exp<-readRDS(hto12_combined_obj_rds)

resolution=c(1:9) / 100
hto.exp <- FindClusters(hto.exp, reduction = "pca", resolution=resolution) 

for(res in resolution){
  colname = paste0("RNA_snn_res.", res)
  g<-DimPlot(hto.exp, reduction = "umap", group.by=colname, label=T) + ggtitle(colname)
  png(paste0("hto12/hto12.exp.res", res, ".png"), width=2000, height=1500, res=300)
  print(g)
  dev.off()
}
