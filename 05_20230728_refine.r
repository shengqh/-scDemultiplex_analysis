source('/home/shengq2/program/scDemultiplex_analysis/common.r')

do_scDemultiplex_refine<-function(root_dir, cur_sample, p.cut=0.001){
  load_install("scDemultiplex", c('shengqh/cutoff', 'shengqh/scDemultiplex'))

  sample_folder=paste0(root_dir, cur_sample)
  setwd(sample_folder)

  result_rds_file = paste0(cur_sample, ".results_obj.rds")
  cat("read", result_rds_file, "\n")
  obj <- readRDS(result_rds_file)

  stopifnot(all(obj$scDemultiplex_full == obj$scDemultiplex))
  
  ntags = nrow(obj)

  result_folder=paste0(sample_folder, "/scDemultiplex_refine")
  if(!dir.exists(result_folder)){
    dir.create(result_folder)
  }
  setwd(result_folder)

  final_file=paste0(cur_sample, ".scDemultiplex_refine.rds")
  
  col=htocols[6]
  for(col in htocols){
    if(grepl("scDemultiplex", col)) {
      next
    }
    cat("refine", col, "\n")

    newcol = paste0(col, "_scDemultiplex")
    prefix = paste0(cur_sample, ".", newcol, "_p", p.cut)
    meta_rds_file = paste0(prefix, ".rds")
    if(file.exists(meta_rds_file)){
      cat("  read ", meta_rds_file, " ...\n")
      meta=readRDS(meta_rds_file)
      obj@meta.data[,newcol]=meta[,newcol]
      obj@meta.data[,paste0(newcol, ".global")]=meta[,paste0(newcol, ".global")]
    }else{
      cat("  scDemultiplex refinement of", col, " ...\n")
      if(!all(rownames(obj) %in% unlist(obj@meta.data[,col]))){
        print("Not all tags in init result, ignored.")
        #missing tags in init result, ignored.
        next;
      }
      obj<-demulti_refine(obj, p.cut, refine_negative_doublet_only=FALSE, mc.cores=ntags, init_column=col)
      obj@meta.data[,newcol]=obj$scDemultiplex
      obj@meta.data[,paste0(newcol, ".global")]=obj$scDemultiplex.global

      saveRDS(obj@meta.data, meta_rds_file)

      obj<-hto_plot(obj, prefix, group.by=newcol)
    }
  }

  saveRDS(obj, final_file)
}

cur_sample = "hto12"
for(cur_sample in samples){
  do_scDemultiplex_refine(root_dir, cur_sample)
}
