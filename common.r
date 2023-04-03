library(scDemultiplex)
library(zoo)
library("R.utils")
library(reshape2)
library(Matrix)
library(data.table)
library(Seurat)
library(tictoc)

root_dir="/nobackup/h_cqs/collaboration/20230301_scrna_hto/"
setwd(root_dir)

samples = c("hto12", "pbmc")

sample_tags = list(
  "hto12" = "HEK_A,HEK_B,HEK_C,K562_A,K562_B,K562_C,KG1_A,KG1_B,KG1_C,THP1_A,THP1_B,THP1_C",
  "pbmc" = "HTO_A,HTO_B,HTO_C,HTO_D,HTO_E,HTO_F,HTO_G,HTO_H"
)

