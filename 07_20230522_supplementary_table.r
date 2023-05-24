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
  # 
  # 
  #   
  #   
  #   for(ftype in c("cell", "ari", "fscore")) {
  # 
  #     stopifnot(all(original_table$HTO == rownames(refine_table)))
  # 
  #     improvement_tbl = cbind(original_table, refine_table)
  # 
  #     if(ftype == "de"){
  #       improvement_tbl = keep_digits(improvement_tbl, 2)
  #     }
  # 
  #     columns = htocols[2:length(htocols)]
  #     col = columns[1]
  #     for(col in columns){
  #       original_col = paste0(col, "_original")
  #       refine_col = paste0(col, "_scDemultiplex")
  #       sign = ifelse(improvement_tbl[,refine_col] > improvement_tbl[,original_col], "↑", "↓")
  #       improvement_tbl[,col] = paste0(improvement_tbl[,refine_col], "(", improvement_tbl[,original_col], ")", sign)
  #     }
  # 
  #     supplementary_tbl = improvement_tbl[,c("HTO", columns)]
  # 
  #     wb <- createWorkbook() # create a workbook
  #     addWorksheet(wb, "Sheet", gridLines = TRUE) #add a worksheet to the workbook
  #     writeData(wb, "Sheet", supplementary_tbl)
  # 
  #     upStyle <- createStyle(fontColour = "red")
  #     downStyle <- createStyle(fontColour = "blue")
  # 
  #     conditionalFormatting(wb, "Sheet", cols = 1:ncol(supplementary_tbl),
  #                   rows = 1:(nrow(supplementary_tbl)+1), rule = "↑", style = upStyle,
  #                   type = "contains")
  # 
  #     conditionalFormatting(wb, "Sheet", cols = 1:ncol(supplementary_tbl),
  #                   rows = 1:(nrow(supplementary_tbl)+1), rule = "↓", style = downStyle,
  #                   type = "contains")
  # 
  #     setColWidths(wb, "Sheet", cols = 1:ncol(supplementary_tbl), widths=c(9, 11, 15, 13, 10, 12) + 2)
  # 
  #     saveWorkbook(wb, paste0("supplementary_table.", name, ".", ftype, ".xlsx"), overwrite = TRUE)
  #   }
  # }
