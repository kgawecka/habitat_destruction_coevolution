# PATCH NETWORK MODULARITY

# - this script computes modularity of local networks
# - for models: v7

library(igraph)
library(data.table)


# function to calculate modularity using IGRAPH package
patch_modularity = function(networkName, networkType, version, n, dD) {
  
  # network incidence matrix
  M_inc = read.table(paste0("../../Data/",networkName,"/M_inc.csv"))
  
  # number of plants, animals and interactions
  n_p = nrow(M_inc)
  n_a = ncol(M_inc)
  n_int = sum(M_inc)
  
  # habitat loss vector
  D = seq(0,1,dD)
  
  # number of patches
  n_patch = n*n
  
  # initialise dataframe for storing results
  df_mod = data.frame(expand.grid(rep_no=1:5, D=D, patch_no=1:n_patch),
                      mod = rep(NA, 5*n_patch*length(D)))  # modularity
  
  for(r in 1:5){
    
    for(k in 1:(length(D)-1)){
      
      # import simulation results
      df_p_plants = fread(paste0("../../Output_",version,"/",networkName,"/",networkType,"_out_",networkName,"_p_plants_r",r,"_D",D[k]*100,".csv"))
      df_p_animals = fread(paste0("../../Output_",version,"/",networkName,"/",networkType,"_out_",networkName,"_p_animals_r",r,"_D",D[k]*100,".csv"))
      
      for(patch in 1:n_patch) {
        
        # data for current patch
        p_current = df_p_plants[df_p_plants$patch_no==patch,"p"]$p
        a_current = df_p_animals[df_p_animals$patch_no==patch,"p"]$p
        n_p_current = sum(p_current)
        n_a_current = sum(a_current)
        
        # if patch contains less than 2 plants or 2 animals, move to next
        if(n_p_current<2 || n_a_current<2) { next }
        
        # current incidence matrix
        M_patch = matrix(1,n_p,n_a) * p_current * rep(a_current, rep(n_p,n_a)) * M_inc
        
        # remove empty rows and columns
        M_patch = M_patch[rowSums(M_patch)>0,colSums(M_patch)>0]
        M_patch = as.matrix(M_patch)
        
        # convert to igraph
        graph = graph_from_incidence_matrix(M_patch)
        
        # determine modules
        modules = cluster_louvain(graph)
        
        # calculate modularity
        Q = modularity(modules)
        
        # store results
        df_mod[df_mod$rep_no==r & df_mod$D==D[k] & df_mod$patch_no==patch,"mod"] = Q
        
      } # patch loop
    } # k loop
  } # r loop
  
  # remove empty rows
  df_mod = df_mod[complete.cases(df_mod), ]
  
  # write out results
  write.csv(df_mod, paste0("../../Results_",version,"/",networkName,"/patch_modularity.csv"), row.names=FALSE)
  
  return(0)
}

# run function
patch_modularity(networkName="M_PL_010", networkType="mutualistic", version="v7", n=100, dD=0.1)
