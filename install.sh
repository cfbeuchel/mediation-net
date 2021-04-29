#!/bin/bash

R -q -e "install.packages(c('visNetwork'), repos = 'https://cran.rstudio.com', quiet = TRUE)"
R -q -e "install.packages(c('magrittr'), quiet = TRUE)"
R -q -e "install.packages(c('data.table'), quiet = TRUE)"

