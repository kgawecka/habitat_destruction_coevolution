# COMBINE STEADY STATE RESULTS

# - this script combines outputs for all values of habitat loss into a single dataframe
# - for models: v9a

# Packages
@everywhere using DataFrames, StatsBase, Random, Distributions, CSV, Dates, DelimitedFiles, LinearAlgebra


@everywhere function combine_results(networkName, networkType, version, n, dD, replicas)

# network incidence matrix
M_inc = readdlm(string("../../Data/",networkName,"/M_inc.csv"), ' ', Int)

# number of plant and animas species
n_p = size(M_inc, 1)
n_a = size(M_inc, 2)

# habitat loss
D = collect(0:dD:1)

# Create empty dataframe
species_presence = vcat(DataFrame(network = fill(networkName, replicas*length(D)*n*n*n_p),
                                  rep_no = repeat(1:replicas, inner=length(D)*n*n*n_p),
                                  D=repeat(D, outer=replicas, inner=(n*n*n_p)),
                                  patch_no=repeat(1:(n*n), outer=(replicas*length(D)*n_p)),
                                  species=repeat(1:n_p, outer=replicas*length(D), inner=(n*n)),
                                  level=fill("resources", replicas*length(D)*n*n*n_p),
                                  p=Vector{Union{Missing, Int64}}(missing, replicas*length(D)*n*n*n_p)),
                        DataFrame(network = fill(networkName, replicas*length(D)*n*n*n_a),
                                  rep_no = repeat(1:replicas, inner=length(D)*n*n*n_a),
                                  D=repeat(D, outer=replicas, inner=(n*n*n_a)),
                                  patch_no=repeat(1:(n*n), outer=(replicas*length(D)*n_a)),
                                  species=repeat(1:n_a, outer=replicas*length(D), inner=(n*n)),
                                  level=fill("consumers", replicas*length(D)*n*n*n_a),
                                  p=Vector{Union{Missing, Int64}}(missing, replicas*length(D)*n*n*n_a)))

for r=1:replicas

for k=1:(length(D)-1)

# Import simulation results
df_p_plants = CSV.read(string("../../Output_",version,"/",networkName,"/",networkType,"_out_",networkName,"_p_plants_r",r,"_D",floor(Int,D[k]*100),".csv"), DataFrame)
df_p_animals = CSV.read(string("../../Output_",version,"/",networkName,"/",networkType,"_out_",networkName,"_p_animals_r",r,"_D",floor(Int,D[k]*100),".csv"), DataFrame)

# Store in dataframe
species_presence[(species_presence.rep_no.===r) .& (species_presence.D.===D[k]) .& (species_presence.level.==="resources"),"p"] .= df_p_plants[:,"p"]
species_presence[(species_presence.rep_no.===r) .& (species_presence.D.===D[k]) .& (species_presence.level.==="consumers"),"p"] .= df_p_animals[:,"p"]

end

end

# Write out results
CSV.write(string("../../Results_",version,"/",networkName,"/p_eq.csv"), species_presence)

return 0
end


# Function for running networks in parallel
@everywhere function run_parallel(network, type, model_version, reps)

  networkName = network
  networkType = type
  version = model_version
  replicas = reps
  n = 100
  dD = 0.1

  combine_results(networkName, networkType, version, n, dD, replicas)

end

# Run simulations
#pmap(run_parallel, [networks], repeat([type], length(networks)), repeat([version], length(networks)))
pmap(run_parallel, ["M_SD_012"], repeat(["mutualistic"], 1), repeat(["v7a"], 1), repeat([5], 1))
#pmap(run_parallel, ["A_PH_004"], repeat(["antagonistic"], 1), repeat(["v7a"], 1), repeat([5], 1))

#################
# TO RUN THIS SCRIPT
# REPLACE NETWORK NAMES, TYPES AND VERSION IN LINE 126
# DO THE FOLLOWING IN THE COMMAND LINE
# julia -p 5 combine_results_replicas_a.jl
