cd /nobackup/h_cqs/collaboration/20230522_scrna_hto

cp -f /home/shengq2/program/scDemultiplex_analysis/03_20230522_check_result.rmd 20230526_check_result.rmd
R -e 'rmarkdown::render("20230526_check_result.rmd")'
