# which packages do I need for this app
# necessary.packages <- c(
#   "shiny", 
#   "magrittr",
#   "data.table",
#   "visNetwork")

# check for installed packages and install them if necessary
# installed <- installed.packages()
# needed <- necessary.packages
# to.install <- needed[!(needed %in% installed[,1])]
# rm(installed)

# if(length(to.install)!=0){
#   if("sva" %in% to.install){
#     install.packages("BiocManager")
# 
# 
#     BiocManager::install("sva")
#   }
#   install.packages(to.install[to.install!="sva"], repos = "https://cran.uni-muenster.de/")
# }

# seperate calls to all packages? -> package dependency is otherwise not recognized by shinyapps.io
library("shiny")
library("data.table")
library("igraph")
library("scales")

# Functions ---------------------------------------------------------------
source("functions/create_mediation_network_data.R")
source("functions/mediation_network.R")
source("functions/allDuplicatedEntries.R")
