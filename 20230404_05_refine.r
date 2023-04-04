source('/home/shengq2/program/scDemultiplex_analysis/common.r')

do_scDemultiplex_refine<-function(root_dir, cur_sample, p.cut=0.001){
  sample_folder=paste0(root_dir, cur_sample)
  setwd(sample_folder)

  result_rds_file = paste0(cur_sample, ".results_obj.rds")
  obj <- readRDS(result_rds_file)

  stopifnot(all(obj$scDemultiplex_full == obj$scDemultiplex))
  
  ntags = nrow(obj)

  result_folder=paste0(sample_folder, "/scDemultiplex_refine")
  if(!dir.exists(result_folder)){
    dir.create(result_folder)
  }
  setwd(result_folder)

  final_file=paste0(cur_sample, ".scDemultiplex_refine.rds")
  
  for(col in htocols){
    if(col == "scDemultiplex") {
      next
    }

    newcol = paste0(col, "_scDemultiplex")
    meta_rds_file = paste0(cur_sample, ".", newcol, "_p", p.cut, ".rds")
    if(file.exists(meta_rds_file)){
      cat("  read ", meta_rds_file, " ...\n")
      obj@meta.data=readRDS(meta_rds_file)
    }else{
      cat("  scDemultiplex refinement of ", col, " ...\n")
      obj<-demulti_refine(obj, p.cut, refine_negative_doublet_only=FALSE, mc.cores=ntags, init_column=col)
      obj@meta.data[,newcol]=obj$scDemultiplex
      obj@meta.data[,paste0(newcol, ".global")]=obj$scDemultiplex.global

      saveRDS(obj@meta.data, paste0(cur_sample, ".", newcol, "_p", p.cut, ".rds"))
    }

    obj<-hto_plot(obj, paste0(cur_sample, ".", newcol, "_p", p.cut), group.by=newcol)
  }

  saveRDS(obj, final_file)
}

cur_sample = "hto12"
for(cur_sample in samples){
  do_scDemultiplex_refine(root_dir, cur_sample)
}
