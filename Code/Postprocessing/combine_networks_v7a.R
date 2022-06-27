# COMBINE ALL NETWORK RESULTS

# - this script combines results of all networks into a single daraframe
# - for models: v7a

# load packages
library(dplyr)
library(ggplot2)
library(tidyr)
library(data.table)
library(FactoMineR)

# set working directory to "habitat_destruction_coevolution"

networksM = c("M_SD_005", "M_SD_008", "M_SD_010", "M_SD_012", "M_SD_025", 
              "M_PL_006", "M_PL_010", "M_PL_036", "M_PL_037", "M_PL_059")

networksA = c("A_HP_008", "A_HP_015", "A_HP_028", "A_HP_032", "A_HP_035", 
              "A_HP_050", "A_PH_004", "A_PH_005", "A_FW_001", "A_FW_006")



# TRAIT MATCHING ----

matching_species_all = data.frame(level=factor(),
                                  species=integer(),
                                  match=double(),
                                  network=factor(),
                                  type=factor(),
                                  version=factor())

for(network in networksM){
  
  # import data
  M_inc = read.table(paste0("Data/",network,"/M_inc.csv"))
  n_p = nrow(M_inc)
  n_a = ncol(M_inc)
  
  M_adj = cbind(rbind(matrix(0, n_p, n_p), t(as.matrix(M_inc))), 
                      rbind(as.matrix(M_inc), matrix(0,n_a,n_a)))
  # Assign NA to 
  M_adj[M_adj==0] = NA
  
  df_z = read.csv(paste0("Data/",network,"/coevolution_z.csv"))
  z = df_z$z
  theta = df_z$theta
  
  alpha = 0.2
  
  # theta matching - v7a
  matching_diff = dist(theta,upper=TRUE,diag=TRUE)
  matching = exp(-alpha*(matching_diff)^2)
  
  # Assign NA to matching between non-interacting species
  matching = as.matrix(matching) * M_adj
  
  # Compute the degree of trait matching for each species
  species_matching = rowMeans(matching, na.rm=TRUE)
  
  # store results
  df_match = data.frame(level=c(rep("resources", n_p), rep("consumers", n_a)),
                        species=c(1:n_p, 1:n_a),
                        match=species_matching,
                        network=as.factor(network),
                        type=as.factor("mutualism"),
                        version="v7a")
  
  matching_species_all = rbind(matching_species_all, df_match)
  
}

for(network in networksA){
  
  # import data
  M_inc = read.table(paste0("Data/",network,"/M_inc.csv"))
  n_p = nrow(M_inc)
  n_a = ncol(M_inc)
  
  M_adj = cbind(rbind(matrix(0, n_p, n_p), t(as.matrix(M_inc))), 
                rbind(as.matrix(M_inc), matrix(0,n_a,n_a)))
  # Assign NA to 
  M_adj[M_adj==0] = NA
  
  df_z = read.csv(paste0("Data/",network,"/coevolution_z.csv"))
  z = df_z$z
  theta = df_z$theta
  
  alpha = 0.2
  
  # theta matching
  matching_diff = dist(theta,upper=TRUE,diag=TRUE)
  matching = exp(-alpha*(matching_diff)^2)
  
  # Assign NA to matching between non-interacting species
  matching = as.matrix(matching) * M_adj
  
  # Compute the degree of trait matching for each species
  species_matching = rowMeans(matching, na.rm=TRUE)
  
  # store results
  df_match = data.frame(level=c(rep("resources", n_p), rep("consumers", n_a)),
                        species=c(1:n_p, 1:n_a),
                        match=species_matching,
                        network=as.factor(network),
                        type=as.factor("antagonism"),
                        version="v7a")
  
  matching_species_all = rbind(matching_species_all, df_match)
  
}

# write out results
write.csv(matching_species_all, paste0("results_combined_v7a/matching_species_all.csv", sep=""), row.names=FALSE)



# ABUNDANCE ----

abundance = data.frame(rep_no=integer(),
                       D=double(),
                       level=factor(),
                       species=factor(),
                       abundance=double(),
                       network=factor(),
                       type=factor(),
                       version=factor())

for(network in networksM){
  
  # import data and calculate abundance
  p_eq = fread(paste0("Results_v7a/",network,"/p_eq.csv")) %>% 
    mutate(species=as.factor(species), level=as.factor(level)) %>% 
    mutate_at("p", ~replace(., is.na(.), 0)) %>%
    group_by(rep_no, D, level, species) %>%
    summarise(abundance = sum(p)/10000) %>%
    ungroup() %>%
    mutate(network=as.factor(network), type=as.factor("mutualism"), version=as.factor("v7a"))
  
  abundance = rbind(abundance, p_eq)
  
}

for(network in networksA){
  
  # import data and calculate abundance
  p_eq = fread(paste0("Results_v7a/",network,"/p_eq.csv")) %>% 
    mutate(species=as.factor(species), level=as.factor(level)) %>% 
    mutate_at("p", ~replace(., is.na(.), 0)) %>%
    group_by(rep_no, D, level, species) %>%
    summarise(abundance = sum(p)/10000) %>%
    ungroup() %>%
    mutate(network=as.factor(network), type=as.factor("antagonism"), version=as.factor("v7a"))
  
  abundance = rbind(abundance, p_eq)
  
}

# write out results
write.csv(abundance, paste0("results_combined_v7a/abundance_v7a.csv", sep=""), row.names=FALSE)



# INTERACTION ABUNDANCE ----

abundance = data.frame(rep_no=integer(),
                       D=double(),
                       int_no=factor(),
                       abundance=double(),
                       network=factor(),
                       type=factor(),
                       version=factor())

for(network in networksM){
  
  # import data
  ab = fread(paste0("Results_v7a/",network,"/abundance_interactions.csv")) %>% 
    mutate(int_no=as.factor(int_no), network=as.factor(network), type=as.factor("mutualism"), version=as.factor("v7a"))
  
  abundance = rbind(abundance, ab)
  
}

for(network in networksA){
  
  # import data
  ab = fread(paste0("Results_v7a/",network,"/abundance_interactions.csv")) %>% 
    mutate(int_no=as.factor(int_no), network=as.factor(network), type=as.factor("antagonism"), version=as.factor("v7a"))
  
  abundance = rbind(abundance, ab)
  
}

# write out results
write.csv(abundance, paste0("results_combined_v7a/abundance_interactions.csv", sep=""), row.names=FALSE)

