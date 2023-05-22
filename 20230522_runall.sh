mkdir -p /nobackup/h_cqs/collaboration/20230522_scrna_hto
cd /nobackup/h_cqs/collaboration/20230522_scrna_hto

if [[ ! -s hashtag-demux-paper ]];then
  git clone https://github.com/Oshlack/hashtag-demux-paper
fi

R -f /home/shengq2/program/scDemultiplex_analysis/00_20230522_benchmark_dataset.r

R -f /home/shengq2/program/scDemultiplex_analysis/01_20230522_analysis.r

R -f /home/shengq2/program/scDemultiplex_analysis/02_20230522_combine_results.r

R -f /home/shengq2/program/scDemultiplex_analysis/03_20230522_check_result.r

cp -f /home/shengq2/program/scDemultiplex_analysis/03_20230522_check_result.rmd 20230522_check_result.rmd
R -e 'rmarkdown::render("20230522_check_result.rmd")'

cp -f /home/shengq2/program/scDemultiplex_analysis/04_20230522_scDemultiplex_iteration.rmd 20230522_scDemultiplex_iteration.rmd
R -e 'rmarkdown::render("20230522_scDemultiplex_iteration.rmd")'

R -f /home/shengq2/program/scDemultiplex_analysis/05_20230522_refine.r

R -f /home/shengq2/program/scDemultiplex_analysis/06_20230522_check_refine_result.r

cp -f /home/shengq2/program/scDemultiplex_analysis/06_20230522_check_refine_result.rmd 20230522_check_refine_result.rmd
R -e 'rmarkdown::render("20230522_check_refine_result.rmd")'

