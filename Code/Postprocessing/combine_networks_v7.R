# COMBINE ALL NETWORK RESULTS

# - this script combines results of all networks into a single daraframe
# - for models: v7

# load packages
library(dplyr)
library(ggplot2)
library(tidyr)
library(data.table)
library(FactoMineR)

# set working directory to "habitat_destruction_coevolution"

networksM = c("M_SD_005", "M_SD_008", "M_SD_010", "M_SD_012", "M_SD_025", 
              "M_PL_006", "M_PL_010", "M_PL_036", "M_PL_037", "M_PL_059")

networksA = c("A_HP_005", "A_HP_008", "A_HP_015", "A_HP_028", "A_HP_032",  
              "A_HP_035", "A_HP_042", "A_HP_050", "A_PH_004", "A_PH_005")


# PATCH NETWORKS ----

patch_networks_all = data.frame(network=factor(),
                                rep_no=integer(),
                                D=double(),
                                patch_no=integer(),
                                state=integer(),
                                int=double(),
                                con=double(),
                                nest=double(),
                                n_sp=integer(),
                                fract_sp=double(),
                                type=factor(),
                                mod=double())

for(network in networksM){
  
  # import data
  M_inc = read.table(paste0("Data/",network,"/M_inc.csv"))
  n_species = nrow(M_inc) + ncol(M_inc)
  
  patch_networks = fread(paste0("Results_v7/",network,"/patch_networks.csv")) %>% 
    filter(n_p>0 & n_a>0) %>%
    mutate(n_sp=(n_p+n_a), fract_sp=(n_p+n_a)/n_species, network=as.factor(network), type=as.factor("mutualism")) %>%
    select(network, rep_no, D, patch_no, state, int, con, nest, n_sp, fract_sp, type)
  
  patch_networks = left_join(patch_networks, fread(paste0(network,"/Results_v7/patch_modularity.csv")))
  
  patch_networks_all = rbind(patch_networks_all, patch_networks)
  
}

for(network in networksA){
  
  # import data
  M_inc = read.table(paste0("Data/",network,"/M_inc.csv"))
  n_species = nrow(M_inc) + ncol(M_inc)
  
  patch_networks = fread(paste0("Results_v7/",network,"/patch_networks.csv")) %>% 
    filter(n_p>0 & n_a>0) %>%
    mutate(n_sp=(n_p+n_a), fract_sp=(n_p+n_a)/n_species, network=as.factor(network), type=as.factor("antagonism")) %>%
    select(network, rep_no, D, patch_no, state, int, con, nest, n_sp, fract_sp, type)
  
  patch_networks = left_join(patch_networks, fread(paste0(network,"/Results_v7/patch_modularity.csv")))
  
  patch_networks_all = rbind(patch_networks_all, patch_networks)
  
}

# write out results
write.csv(patch_networks_all, paste0("results_combined_v7/patch_networks_all.csv", sep=""), row.names=FALSE)



# PATCH NEIGHBOURS ----

patch_neighbours_all = data.frame(network=factor(),
                                  rep_no=integer(),
                                  D=double(),
                                  patch_no=integer(),
                                  no_neigh=integer(),
                                  no_occ_neigh=integer(),
                                  no_species=integer(),
                                  type=factor())

for(network in networksM){
  
  # import data
  patch_neigh = fread(paste0("Results_v7/",network,"/patch_neighbours.csv")) %>% 
    drop_na %>%
    mutate(network=as.factor(network), type=as.factor("mutualism"))
  
  patch_neighbours_all = rbind(patch_neighbours_all, patch_neigh)
  
}

for(network in networksA){
  
  # import data
  patch_neigh = fread(paste0("Results_v7/",network,"/patch_neighbours.csv")) %>% 
    drop_na %>%
    mutate(network=as.factor(network), type=as.factor("antagonism"))
  
  patch_neighbours_all = rbind(patch_neighbours_all, patch_neigh)
  
}

# write out results
write.csv(patch_neighbours_all, paste0("results_combined_v7/patch_neighbours_all.csv", sep=""), row.names=FALSE)



# NETWORK STRUCTURE PCA ----

# import data
patch_networks_all = fread(paste0("results_combined_v7/patch_networks_all.csv"))

# remove NAs
patch_networks_pca = patch_networks_all[complete.cases(patch_networks_all), ]

# PCA
pca_fsp_con_nest_mod = PCA(patch_networks_pca %>% select(fract_sp, con, nest, mod), graph=TRUE)
pca_nsp_con_nest_mod = PCA(patch_networks_pca %>% select(n_sp, con, nest, mod), graph=TRUE)

# store results
patch_networks_pca_all = patch_networks_pca %>%
  select(type, network, rep_no, D, patch_no) %>%
  mutate(PC1_fsp_con_nest_mod=as.vector(pca_fsp_con_nest_mod$ind$coord[,1]), 
         PC2_fsp_con_nest_mod=as.vector(pca_fsp_con_nest_mod$ind$coord[,2]),
         PC1_nsp_con_nest_mod=as.vector(pca_nsp_con_nest_mod$ind$coord[,1]), 
         PC2_nsp_con_nest_mod=as.vector(pca_nsp_con_nest_mod$ind$coord[,2]))

# write out results
write.csv(patch_networks_pca_all, paste0("results_combined_v7/patch_networks_pca_all.csv", sep=""), row.names=FALSE)



# TRAIT VALUES ----

z_species = data.frame(D=double(),
                       level=factor(),
                       species=integer(),
                       z_mean=double(),
                       z_median=double(),
                       z_sd=double(),
                       z_IQR=double(),
                       network=factor(),
                       type=factor())

for(network in networksM){
  
  # import data and summarise across replicas, patches and species
  z_summ = fread(paste0("Results_v7/",network,"/z_eq.csv")) %>% 
    group_by(D, level, species) %>%
    summarise(z_mean = mean(z, na.rm=TRUE),
              z_median = median(z, na.rm=TRUE),
              z_sd = sd(z, na.rm=TRUE),
              z_IQR = IQR(z, na.rm=TRUE)) %>%
    ungroup() %>%
    mutate(network=as.factor(network), type=as.factor("mutualism"))
  
  z_species = rbind(z_species, z_summ)
  
}

for(network in networksA){
  
  # import data and summarise across replicas, patches and species
  z_summ = fread(paste0("Results_v7/",network,"/z_eq.csv")) %>% 
    group_by(D, level, species) %>%
    summarise(z_mean = mean(z, na.rm=TRUE),
              z_median = median(z, na.rm=TRUE),
              z_sd = sd(z, na.rm=TRUE),
              z_IQR = IQR(z, na.rm=TRUE)) %>%
    ungroup() %>%
    mutate(network=as.factor(network), type=as.factor("antagonism"))
  
  z_species = rbind(z_species, z_summ)
  
}

# write out results
write.csv(z_species, paste0("results_combined_v7/z_species.csv", sep=""), row.names=FALSE)



# TRAIT MATCHING ----

matching_network_all = data.frame(network=factor(),
                                  rep_no=integer(),
                                  D=double(),
                                  patch_no=integer(),
                                  match=double(),
                                  match_theta=double(),
                                  match_resources=double(),
                                  match_consumers=double(),
                                  type=factor())

for(network in networksM){
  
  matching = fread(paste0("Results_v7/",network,"/matching_network.csv")) %>% 
    drop_na %>%
    mutate(network=as.factor(network), type=as.factor("mutualism"))
  
  matching_network_all = rbind(matching_network_all, matching)
  
}

for(network in networksA){
  
  matching = fread(paste0("Results_v7/",network,"/matching_network.csv")) %>% 
    drop_na %>%
    mutate(network=as.factor(network), type=as.factor("antagonism"))
  
  matching_network_all = rbind(matching_network_all, matching)
  
}

# write out results
write.csv(matching_network_all, paste0("results_combined_v7/matching_network_all.csv", sep=""), row.names=FALSE)



# INTERACTION ABUNDANCE ----

abundance = data.frame(rep_no=integer(),
                       D=double(),
                       int_no=factor(),
                       abundance=double(),
                       network=factor(),
                       type=factor())

for(network in networksM){
  
  # import data
  ab = fread(paste0("Results_v7/",network,"/abundance_interactions.csv")) %>% 
    mutate(int_no=as.factor(int_no), network=as.factor(network), type=as.factor("mutualism"))
  
  abundance = rbind(abundance, ab)
  
}

for(network in networksA){
  
  # import data
  ab = fread(paste0("Results_v7/",network,"/abundance_interactions.csv")) %>% 
    mutate(int_no=as.factor(int_no), network=as.factor(network), type=as.factor("antagonism"))
  
  abundance = rbind(abundance, ab)
  
}

# write out results
write.csv(abundance, paste0("results_combined_v7/abundance_interactions.csv", sep=""), row.names=FALSE)



# ABUNDANCE ----

abundance = data.frame(rep_no=integer(),
                       D=double(),
                       level=factor(),
                       species=factor(),
                       abundance=double(),
                       network=factor(),
                       type=factor())

for(network in networksM){
  
  # import data and calculate abundance
  p_eq = fread(paste0("Results_v7/",network,"/p_eq.csv")) %>% 
    mutate(species=as.factor(species), level=as.factor(level)) %>% 
    mutate_at("p", ~replace(., is.na(.), 0)) %>%
    group_by(rep_no, D, level, species) %>%
    summarise(abundance = sum(p)/10000) %>%
    ungroup() %>%
    mutate(network=as.factor(network), type=as.factor("mutualism"))
  
  abundance = rbind(abundance, p_eq)
  
}

for(network in networksA){
  
  # import data and calculate abundance
  p_eq = fread(paste0("Results_v7/",network,"/p_eq.csv")) %>% 
    mutate(species=as.factor(species), level=as.factor(level)) %>% 
    mutate_at("p", ~replace(., is.na(.), 0)) %>%
    group_by(rep_no, D, level, species) %>%
    summarise(abundance = sum(p)/10000) %>%
    ungroup() %>%
    mutate(network=as.factor(network), type=as.factor("antagonism"))
  
  abundance = rbind(abundance, p_eq)
  
}

# write out results
write.csv(abundance, paste0("results_combined_v7/abundance.csv", sep=""), row.names=FALSE)


# SPECIES DEGREE ----

species_degree = data.frame(type=factor(), 
                            network=factor(),
                            level=factor(),
                            species=integer(),
                            degree=integer(),
                            degree_fract=double())

for(network in networksM){
  
  # import data
  M_inc = read.table(paste0("Data/",network,"/M_inc.csv"))
  n_p = nrow(M_inc)
  n_a = ncol(M_inc)
  n_sp = n_p + n_a
  
  df_degree = data.frame(type=rep("mutualism", n_sp),
                         network=rep(network, n_sp),
                         level=c(rep("resources", n_p), rep("consumers", n_a)),
                         species=c(1:n_p, 1:n_a),
                         degree=c(rowSums(M_inc), colSums(M_inc)),
                         degree_fract=c(rowSums(M_inc)/n_a, colSums(M_inc)/n_p))
  
  species_degree = rbind(species_degree, df_degree)
  
}

for(network in networksA){
  
  # import data
  M_inc = read.table(paste0("Data/",network,"/M_inc.csv"))
  n_p = nrow(M_inc)
  n_a = ncol(M_inc)
  n_sp = n_p + n_a
  
  df_degree = data.frame(type=rep("antagonism", n_sp),
                         network=rep(network, n_sp),
                         level=c(rep("resources", n_p), rep("consumers", n_a)),
                         species=c(1:n_p, 1:n_a),
                         degree=c(rowSums(M_inc), colSums(M_inc)),
                         degree_fract=c(rowSums(M_inc)/n_a, colSums(M_inc)/n_p))
  
  species_degree = rbind(species_degree, df_degree)
  
}

# write out results
write.csv(species_degree, paste0("results_combined_v7/species_degree.csv", sep=""), row.names=FALSE)



# PATCH DISSIMILARITY ----

patch_beta_all = data.frame(rep_no=integer(),
                            D=double(),
                            patch_i=integer(),
                            patch_j=integer(),
                            beta_S=double(),
                            beta_OS=double(),
                            beta_WN=double(),
                            beta_ST=double(),
                            network=factor(),
                            type=factor())

for(network in networksM){
  
  # import data
  patch_beta = fread(paste0("Results_v7/",network,"/patch_dissimilarity.csv")) %>% 
    select(rep_no, D, patch_i, patch_j, beta_S, beta_OS, beta_WN, beta_ST) %>%
    mutate(network=as.factor(network), type=as.factor("mutualism"))
  
  patch_beta_all = rbind(patch_beta_all, patch_beta)
  
}

for(network in networksA){
  
  # import data
  patch_beta = fread(paste0("Results_v7/",network,"/patch_dissimilarity.csv")) %>% 
    select(rep_no, D, patch_i, patch_j, beta_S, beta_OS, beta_WN, beta_ST) %>%
    mutate(network=as.factor(network), type=as.factor("antagonism"))
  
  patch_beta_all = rbind(patch_beta_all, patch_beta)
  
}

# write out results
write.csv(patch_beta_all, paste0("results_combined_v7/patch_beta_all.csv", sep=""), row.names=FALSE)
