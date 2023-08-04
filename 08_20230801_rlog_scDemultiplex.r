#install.packages('/home/shengq2/program/scDemultiplex', repos = NULL, type="source")

source('/home/shengq2/program/scDemultiplex_analysis/common.r')

cur_sample = "barnyard"
for(cur_sample in samples){
  do_scDemultiplex(root_dir, cur_sample, p.cuts=scDemultiplex.p.cuts, do_rlog=TRUE)
}
