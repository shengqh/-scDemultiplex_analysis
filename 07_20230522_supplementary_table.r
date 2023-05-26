library(kableExtra)
library(data.table)
library(openxlsx)

source("/home/shengq2/program/scDemultiplex_analysis/common.r")

setwd(root_dir)

keep_digits<-function(DF, digits){
  is.num <- sapply(DF, is.numeric)
  DF[is.num] <- lapply(DF[is.num], round, digits)
  return(DF)
}

cols = htocols[htocols != "scDemultiplex"]


m = "ARI"
for(m in c("ARI", "Fscore")){
  df = data.frame(matrix(NA, ncol=length(samples), nrow = length(cols)))
  colnames(df) <- samples
  rownames(df) <- cols
  
  name=samples[1]
  for(name in samples){
    refine_table = data.frame(fread(file.path(root_dir, name, "scDemultiplex_refine", paste0(name, ".score.csv"))), row.names=1)
    rt = refine_table[rownames(refine_table) != "scDemultiplex",m, drop=F]
    rt = keep_digits(rt, 3)
    
    col = cols[1]
    for (col in cols){
      old_v = rt[col, m]
      new_v = rt[paste0(col, "_scDemultiplex"), m]
      if(is.na(old_v)){
        df[col, name] = ""
      }else{
        df[col, name] = paste0(new_v, " (", old_v, ")")
      }
    }
  }
  write.csv(df, paste0("supplementary_table_", m, ".csv"))
}
