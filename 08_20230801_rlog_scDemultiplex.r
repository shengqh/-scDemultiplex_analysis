#install.packages('/home/shengq2/program/scDemultiplex', repos = NULL, type="source")

source('/home/shengq2/program/scDemultiplex_analysis/common.r')

cur_sample = "barnyard"

#valid_samples = samples[samples != "barnyard"]
valid_samples = samples
for(cur_sample in valid_samples){
  do_scDemultiplex(root_dir, cur_sample, p.cuts=scDemultiplex.p.cuts, do_rlog=TRUE)
}
