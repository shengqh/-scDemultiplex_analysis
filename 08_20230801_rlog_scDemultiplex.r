source('/home/shengq2/program/scDemultiplex_analysis/common.r')

cur_sample = "barnyard"
do_scDemultiplex_rlog<-function(root_dir, cur_sample, p.cut=0.001){
  load_install("scDemultiplex", c('shengqh/cutoff', 'shengqh/scDemultiplex'))

  sample_folder=paste0(root_dir, cur_sample)
  setwd(sample_folder)

  result_rds_file = paste0(cur_sample, ".results_obj.rds")
  cat("read", result_rds_file, "\n")
  obj <- readRDS(result_rds_file)
  obj <- NormalizeData(obj, normalization.method = "LogNormalize")

  ntags = nrow(obj)

  result_folder=paste0(sample_folder, "/scDemultiplex_rlog")
  if(!dir.exists(result_folder)){
    dir.create(result_folder)
  }
  setwd(result_folder)

  output_prefix<-paste0(cur_sample, ".HTO")

  final_file=paste0(cur_sample, ".scDemultiplex_rlog.cutoff.rds")
  if(file.exists(final_file)){
    cat("  ", final_file, " already exists, skip\n")
  }else{
    tic(paste0("starting ", cur_sample, " cutoff using log normalized data...\n"))
    cat("  scDemultiplex_cutoff ...\n")
    obj<-demulti_cutoff(obj, output_prefix=output_prefix, cutoff_startval = 0, mc.cores=ntags)
    toc1=toc()

    cutoffs = readRDS(paste0(output_prefix, ".cutoff_list.rds"))
    failed = names(cutoffs)[is.na(as.numeric(cutoffs))]
    if(length(failed) > 0){
      cat("  failed cutoffs: ", paste0(failed, collapse=", "), "\n")

      failedTag = failed[1]
      for(failedTag in failed){
        values = FetchData(obj, vars = failedTag)[,1]
        png(paste0(output_prefix, "_", failedTag, ".cutoff.png"), width = 2000, height = 1600, res = 300)
        hist(values, 200, F, xlab = "concentration", ylab = "density", main = NULL, col = "grey")
        lines(density(values), lwd = 1.5, col = "blue")
        dev.off()
      }
      writeLines(failed, paste0(output_prefix, ".failed_tags.txt"))
    }else{
      obj<-hto_plot(obj, paste0(output_prefix, ".cutoff"), group.by="scDemultiplex_cutoff")
      saveRDS(obj, final_file)
    }
  }
}

cur_sample = "pbmc8"
for(cur_sample in samples){
  do_scDemultiplex_rlog(root_dir, cur_sample)
}
