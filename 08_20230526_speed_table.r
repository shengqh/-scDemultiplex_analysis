library(kableExtra)
library(data.table)
library(openxlsx)

is_unix=.Platform$OS.type == "unix"
if(is_unix) {
  source('/home/shengq2/program/scDemultiplex_analysis/common.r')
} else {
  source('C:/Users/sheng/Programs/scDemultiplex_analysis/common.r')
}

setwd(root_dir)

cells=unlist(lapply(samples, function(x){
  read.csv(paste0(x, "/", x, ".ncell.csv"))[1,"Cell"]
}))

speed = data.frame(t(read.csv("speed.csv", row.names=1)))

speed$Cells=cells
speed=speed[,c("Cells", "scDemultiplex_cutoff", "scDemultiplex")]

speed$Computer="ACCRE Cluster Gateway"
speed$CPU="Intel(R) Xeon(R) Gold 6154 CPU @ 3.00GHz"
speed$Memory="1000gb"
speed$System="CentOs 7"
speed$R="4.3.0"
speed=speed[,c("Computer", "CPU", "Memory", "System", "R", "Cells", "scDemultiplex_cutoff", "scDemultiplex")]
speed$scDemultiplex_cutoff = paste0(round(speed$scDemultiplex_cutoff, 1), " sec")
speed$scDemultiplex = paste0(round(speed$scDemultiplex/60, 1), " min")

colnames(speed)[7:8] = c("scDemultiplex cutoff", "scDemultiplex refinement")
write.csv(speed, paste0(.Platform$OS.type, ".speed.csv"))
