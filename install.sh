#!/bin/bash

R -q -e "install.packages(c('igraph'), repos = 'https://cran.rstudio.com', quiet = TRUE)"
R -q -e "install.packages(c('data.table'), quiet = TRUE)"
R -q -e "install.packages(c('scales'), quiet = TRUE)"

