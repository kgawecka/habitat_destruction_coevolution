# DOWNLOAD NETWORK

# this script downloads networks from www.web-of-life.es and saves the incidence matrix

library(knitr)
library(bipartite)
library(dplyr)
library(ggplot2)
library(tidyr)
library(bipartite)

# set working directory to "habitat_destruction_coevolution"


# specify network to download
networkName <- "A_HP_050" 


# download network from web of life
url <- paste("http://www.web-of-life.es/download/",
             networkName,
             "_",
             "yes",     # species name?
             ".csv",
             sep = "")
M_inc <- as.matrix(read.csv(file = url, sep = ",", header = TRUE, row.names = 1, check.names = "FALSE"))  # adjacency matrix
M_inc[M_inc > 0] <- 1                      # from quantitative to binary matrix

if("Num. of hosts sampled"%in%colnames(M_inc)){    # remove column "number of hosts sampled"
  denum  <- which(colnames(M_inc)=="Num. of hosts sampled")
  M_inc <- M_inc[,-denum]
}

M_inc <- sortweb(M_inc, sort.order="dec")  # sort network by decreasing number of interactions per species

n_p <- nrow(M_inc)    # no of plants
n_a <- ncol(M_inc)    # no of animals
n_int <- sum(M_inc)   # no of interactions

# write out
write.table(M_inc, paste0("Data/",networkName,"/M_inc.csv"), row.names=FALSE, col.names=FALSE)
