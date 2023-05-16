cd /nobackup/h_cqs/collaboration/20230301_scrna_hto

R -f /home/shengq2/program/scDemultiplex_analysis/20230304_03_check_result.r

cp -f /home/shengq2/program/scDemultiplex_analysis/20230304_04_check_result.rmd 20230514_04_check_result.rmd

#Failed to generate report using R4.3.0. So R4.1.3 was used
R41 -e 'rmarkdown::render("20230514_04_check_result.rmd")'
