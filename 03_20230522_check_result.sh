cd /nobackup/h_cqs/collaboration/20230522_scrna_hto

#R -f /home/shengq2/program/scDemultiplex_analysis/03_20230519_check_result.r

cp -f /home/shengq2/program/scDemultiplex_analysis/03_20230519_check_result.rmd 03_20230519_check_result.rmd

R -e 'rmarkdown::render("03_20230519_check_result.rmd")'
