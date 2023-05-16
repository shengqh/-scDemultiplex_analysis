library(kableExtra)

source("/home/shengq2/program/scDemultiplex_analysis/common.r")

setwd(root_dir)
file.copy("/home/shengq2/program/scDemultiplex_analysis/20230404_07_refine_report.rmd", getwd(), overwrite=TRUE)
rmarkdown::render("20230513_07_refine_report.rmd")