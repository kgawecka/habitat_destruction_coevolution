# REMOVE SPECIES AND CALCULATE NESTEDNESS AND MODULARITY

# species are removed from networks in the order in which they become extinct from the landscape:
#  - MUTUALISM: most specialist species first
#  - ANTAGONISM: most specialist consumers and most generalist resources first
# every time a species is removed, network nestedness and modularity are recalculated

library(igraph)
library(bipartite)
library(dplyr)
library(ggplot2)
library(RColorBrewer)
library(ggpubr)

# nestedness function
nestedness_s <- function(x){ # x is the adjacency matrix (binary)
  # nestedness of rows 
  nested_rows <- 0
  for(i in 1:nrow(x)){
    for(j in i:nrow(x)){
      if(j>i){
        shared <- sum(x[i,]*x[j,]) # sum of common interactions
        min_shared <- min(sum(x[i,]),sum(x[j,])) # min of the degrees
        nested_rows <- nested_rows+(shared/min_shared)
      }
    }
  }
  
  nestedness_rows <- nested_rows/(nrow(x)*(nrow(x)-1)/2)
  
  # nestedness of columns
  nested_columns <- 0
  for(i in 1:ncol(x)){
    for(j in i:ncol(x)){
      if(j>i){
        shared <- sum(x[,i]*x[,j]) # sum of common interactions
        min_shared <- min(sum(x[,i]),sum(x[,j])) # min of the degrees
        nested_columns <- nested_columns+(shared/min_shared)
      }
    }
  }
  nestedness_columns <- nested_columns/(ncol(x)*(ncol(x)-1)/2)
  
  # nestedness of the network
  
  nestedness_network <- (nested_rows+nested_columns)/((nrow(x)*(nrow(x)-1)/2)+(ncol(x)*(ncol(x)-1)/2))
  
  return(nestedness_network)
  
}

# modularity function
measure_modularity <- function(x){
  
  graph = graph_from_incidence_matrix(network)
  modules = cluster_louvain(graph)
  Q = modularity(modules)
  return(Q)
}

network_cleanup <- function(x){
  x[x>0]<-1
  x <- bipartite::empty(x)
}

# remove species function
network_remove_species <- function(x, interaction_type){
  
  if (interaction_type == 'mutualistic') {
    
    row_sum_min_value <- min(rowSums(network))
    col_sum_min_value <- min(colSums(network))
    
    row_sum_min_species <- which(rowSums(network) == row_sum_min_value)
    col_sum_min_species <- which(colSums(network) == col_sum_min_value)
    
    row_sp_to_remove <- sample(x = names(row_sum_min_species), size = 1)
    col_sp_to_remove <- sample(x = names(col_sum_min_species), size = 1)
  }

  if (interaction_type == 'antagonistic'){
    
    row_sum_max_value <- max(rowSums(network))
    col_sum_min_value <- min(colSums(network))
    
    row_sum_max_species <- which(rowSums(network) == row_sum_max_value)
    col_sum_min_species <- which(colSums(network) == col_sum_min_value)
    
    row_sp_to_remove <- sample(x = names(row_sum_max_species), size = 1)
    col_sp_to_remove <- sample(x = names(col_sum_min_species), size = 1)
    
  }
  
  network <- network[rownames(network) != row_sp_to_remove, colnames(network) != col_sp_to_remove]
  
  return(network)
  
}
  
network_size <- function(x){
  size <- nrow(x) + ncol(x)
  return(size)
}

clean_df <- function(x){
  df <- as.data.frame(x)
  rownames(df) <- NULL
  colnames(df) <- c('network_name', 'interaction_type','size','fraction_size','nestedness','modularity')
  return(df)
}


# Mutualistic networks
list_of_mutualistic_networks <- c('M_PL_006','M_PL_010','M_PL_036','M_PL_037','M_PL_059','M_SD_005','M_SD_008','M_SD_010','M_SD_012','M_SD_025')
type <- 'mutualistic'
df_mutualistic <- NULL

for (file in list_of_mutualistic_networks) {
  network <- read.table(paste0("../Data/",file,"/M_inc.csv"))
  name <- tools::file_path_sans_ext(file)
  og_size <- network_size(network)
  nrow_network <- nrow(network)
  ncol_network <- ncol(network)
  print(name)
  
  while((nrow_network != 0) && (ncol_network != 0)){
    network <- network_cleanup(network)
    nestedness <- nestedness_s(network)
    modularity <- measure_modularity(network)
    size <- network_size(network)
    frac_species <- size/og_size
    fill_vector <- c(name, type, size, frac_species, nestedness, modularity)
    df_mutualistic <- rbind(df_mutualistic, fill_vector)
    network <- network_remove_species(network, type)
    if (is.null(dim(network))) {
      nrow_network <- 0
    }
    else{
      nrow_network <- nrow(network)
      ncol_network <- ncol(network)
    }
    
  }

}

df_mutualistic <- clean_df(df_mutualistic)


# antagonistic interactions
list_of_antagonistic_networks <- c('A_HP_005','A_HP_008','A_HP_015','A_HP_028','A_HP_032','A_HP_035','A_HP_042','A_HP_050','A_PH_004','A_PH_005')
type <- 'antagonistic'
df_antagonistic <- NULL

for (file in list_of_antagonistic_networks) {
  network <- read.table(paste0("../Data/",file,"/M_inc.csv"))
  name <- tools::file_path_sans_ext(file)
  og_size <- network_size(network)
  nrow_network <- nrow(network)
  ncol_network <- ncol(network)
  print(name)
  
  while((nrow_network != 0) && (ncol_network != 0)){
    network <- network_cleanup(network)
    nestedness <- nestedness_s(network)
    modularity <- measure_modularity(network)
    size <- network_size(network)
    frac_species <- size/og_size
    fill_vector <- c(name, type, size, frac_species, nestedness, modularity)
    df_antagonistic <- rbind(df_antagonistic, fill_vector)
    network <- network_remove_species(network, type)
    if (is.null(dim(network))) {
      nrow_network <- 0
    }
    else{
      nrow_network <- nrow(network)
      ncol_network <- ncol(network)
    }
    
  }
  
}

# clean dataframe
df_antagonistic <- clean_df(df_antagonistic)

# merge results from antagonism and mutualism
full_df <- rbind(df_mutualistic, df_antagonistic) %>%
  mutate(size=as.numeric(levels(size))[size],
         fraction_size=as.numeric(levels(fraction_size))[fraction_size],
         nestedness=as.numeric(levels(nestedness))[nestedness],
         modularity=as.numeric(levels(modularity))[modularity])

# rename interaction types
full_df$interaction_type2 = factor(full_df$interaction_type, labels=c("mutualism", "antagonism"))

# color palettes with more than 11 colors
getPalette = colorRampPalette(brewer.pal(4, "RdYlBu"))

a = ggplot(data=full_df, aes(x=1-fraction_size, y=nestedness, col=as.factor(network_name))) +
  geom_line() +
  geom_point(alpha=0.3) + 
  facet_grid(~factor(interaction_type2, levels=c("mutualism", "antagonism"))) +
  lims(x=c(0,1), y=c(0,1)) +
  xlab(NULL) + 
  ylab('nestedness') +
  #scale_colour_brewer(palette="RdYlBu") +
  scale_color_manual(values=getPalette(20)) +
  #scale_colour_grey() +
  theme(panel.grid=element_blank(), panel.background=element_blank(), panel.border=element_blank(),
        axis.line=element_line(colour="black"), axis.text=element_text(color="black"), 
        axis.title.x=element_text(vjust=-0.5), axis.title.y=element_text(vjust=2),
        strip.background=element_rect(fill=NA), legend.key=element_blank(), legend.position='none')

b = ggplot(data=full_df, aes(x=1-fraction_size, y=modularity, col=as.factor(network_name))) +
  geom_line() +
  geom_point(alpha=0.3) + 
  facet_grid(~factor(interaction_type2, levels=c("mutualism", "antagonism"))) +
  lims(x=c(0,1), y=c(0,1)) +
  xlab('fraction of species removed') + 
  ylab('modularity') +
  #scale_colour_brewer(palette="RdYlBu") +
  scale_color_manual(values=getPalette(20)) +
  #scale_colour_grey() +
  theme(panel.grid=element_blank(), panel.background=element_blank(), panel.border=element_blank(),
        axis.line=element_line(colour="black"), axis.text=element_text(color="black"), 
        axis.title.x=element_text(vjust=-0.5), axis.title.y=element_text(vjust=2),
        strip.background=element_rect(fill=NA), legend.key=element_blank(), legend.position='none',
        strip.text.x=element_blank())

ggarrange(a, b, ncol=1, nrow=2)
#ggsave("species_removal.pdf", width=160, height=160, units="mm")

