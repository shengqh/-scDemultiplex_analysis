cd /nobackup/h_cqs/collaboration/20230522_scrna_hto

#R -f /home/shengq2/program/scDemultiplex_analysis/06_20230521_check_refine_result.r

cp -f /home/shengq2/program/scDemultiplex_analysis/06_20230521_check_refine_result.rmd 20230522_check_refine_result.rmd

R -e 'rmarkdown::render("20230522_check_refine_result.rmd")'
