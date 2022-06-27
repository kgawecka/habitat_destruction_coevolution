# DISSIMILARITY BETWEEN PATCHES (BETA DIVERSITY)

# - this script computes beta diversity between pairs of patches
# - for models: v7

library(data.table)
library(dplyr)
library(igraph)
library(plyr)
library(reshape2)
library(betalink)
library(bipartite)

# function to calculate beta divversity using BETALINK package
patch_dissimilarity = function(networkName, version, n_sample) {
  
  # network incidence matrix
  M_inc = read.table(paste0("../../Data/",networkName,"/M_inc.csv"))
  M_inc = as.matrix(M_inc)
  
  # number of plants & animals
  n_p = nrow(M_inc)
  n_a = ncol(M_inc)
  
  # import presence data
  p_eq = fread(paste0("../../Results_",version,"/",networkName,"/p_eq.csv"))
  
  # habitat loss vector
  DD = unique(p_eq$D)
  
  # initialise dataframe for storing results
  betadiv = data.frame(expand.grid(rep_no=1:5, D=DD, N=1:(n_sample*(n_sample-1)/2)),
                       patch_i=NA,
                       patch_j=NA,
                       beta_S=NA,
                       beta_OS=NA,
                       beta_WN=NA,
                       beta_ST=NA)
  
  for(r in 1:5){
    
    for(k in 1:(length(DD)-1)){
      
      # filter presnece data for current r and k
      p_D = filter(p_eq, rep_no==r, D==DD[k])
      
      # separate presence data into resources and consumers
      p_res = p_D %>% filter(level=="resources") %>% select(patch_no, species, p) %>% 
        dplyr::rename(species_res=species, p_res=p)
      p_con = p_D %>% filter(level=="consumers") %>% select(patch_no, species, p) %>% 
        dplyr::rename(species_con=species, p_con=p)
      
      # dataframe with all possible interactions
      df = data.frame(expand.grid(species_res=1:n_p, species_con=1:n_a, patch_no=1:10000),
                      M_inc=rep(as.vector(M_inc), 10000))
      
      # filter for possible interactions only, then join with presence data, then filter for species present
      df_int = filter(df, M_inc==1) %>% 
        left_join(., p_res, by=c("patch_no", "species_res")) %>% filter(., p_res==1) %>%
        left_join(., p_con, by=c("patch_no", "species_con")) %>% filter(., p_con==1)
      
      if(nrow(df_int)==0) { next }
      
      patches = unique(df_int$patch_no)
      
      if(length(patches)<=1) { next }
      
      # convert species numbers to factors (must be as factors for betadiversity)
      df_int = df_int %>%
        mutate(species_res = as.factor(paste0("r",species_res)),
               species_con = as.factor(paste0("c",species_con)))
      
      if(length(patches)>n_sample){
        # sample patches
        p_samp = sample(patches, n_sample)
        df_int_D = df_int %>% filter(patch_no %in% p_samp)
      } else {
        df_int_D = df_int
      }
      
      # create list of local networks
      nets = dlply(df_int_D, .(patch_no), acast, species_res~species_con, sum, value.var="M_inc")
      networks = prepare_networks(nets, directed=TRUE)
      
      # compute beta diversity
      betadiv_D = network_betadiversity(networks)

      no_res = nrow(betadiv_D)
      
      # store results
      betadiv[betadiv$rep_no==r & betadiv$D==DD[k],"patch_i"][1:no_res] = as.integer(levels(betadiv_D$i))[betadiv_D$i]
      betadiv[betadiv$rep_no==r & betadiv$D==DD[k],"patch_j"][1:no_res] = as.integer(levels(betadiv_D$j))[betadiv_D$j]
      betadiv[betadiv$rep_no==r & betadiv$D==DD[k],"beta_S"][1:no_res] = betadiv_D$S
      betadiv[betadiv$rep_no==r & betadiv$D==DD[k],"beta_OS"][1:no_res] = betadiv_D$OS
      betadiv[betadiv$rep_no==r & betadiv$D==DD[k],"beta_WN"][1:no_res] = betadiv_D$WN
      betadiv[betadiv$rep_no==r & betadiv$D==DD[k],"beta_ST"][1:no_res] = betadiv_D$ST
      
    } # k loop
  } # r loop
  
  # remove empty rows
  betadiv=filter(betadiv, !is.na(patch_i))
  
  # write out results
  write.csv(betadiv, paste0("../../Results_",version,"/",networkName,"/patch_dissimilarity.csv"), row.names=FALSE)
  
  return(0)
  
}

# run function
patch_dissimilarity(networkName="A_PH_004", version="v7", n_sample=50)
