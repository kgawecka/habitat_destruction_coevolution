# COMBINE STEADY STATE RESULTS

# - this script combines outputs for all values of habitat loss into a single dataframe
# - for models: v7

# Packages
@everywhere using DataFrames, StatsBase, Random, Distributions, CSV, Dates, DelimitedFiles, LinearAlgebra


@everywhere function combine_results(networkName, networkType, version, n, dD, replicas)

# network incidence matrix
M_inc = readdlm(string("../../Data/",networkName,"/M_inc.csv"), ',', Int)

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

z_values = DataFrame(rep_no = Int64[],
                     D = Float64[],
                     patch_no = Int64[],
                     species = Int64[],
                     level = String[],
                     z = Float64[])

ab_traits_dt = vcat(DataFrame(rep_no = repeat(1:replicas, inner=length(D)*1000*n_p),
                       D = repeat(D, outer=replicas, inner=1000*n_p),
                       t = repeat(1:1000, outer=replicas*length(D), inner=n_p),
                       species=repeat(1:n_p, outer=replicas*length(D)*1000),
                       level=fill("resources", replicas*length(D)*1000*n_p),
                       ab=Vector{Union{Missing, Float32}}(missing, replicas*length(D)*1000*n_p),
                       z_mean=Vector{Union{Missing, Float32}}(missing, replicas*length(D)*1000*n_p),
                       z_sd=Vector{Union{Missing, Float32}}(missing, replicas*length(D)*1000*n_p)),
                    DataFrame(rep_no = repeat(1:replicas, inner=length(D)*1000*n_a),
                       D = repeat(D, outer=replicas, inner=1000*n_a),
                       t = repeat(1:1000, outer=replicas*length(D), inner=n_a),
                       species=repeat(1:n_a, outer=replicas*length(D)*1000),
                       level=fill("consumers", replicas*length(D)*1000*n_a),
                       ab=Vector{Union{Missing, Float32}}(missing, replicas*length(D)*1000*n_a),
                       z_mean=Vector{Union{Missing, Float32}}(missing, replicas*length(D)*1000*n_a),
                       z_sd=Vector{Union{Missing, Float32}}(missing, replicas*length(D)*1000*n_a)))

patch_state = DataFrame(rep_no = repeat(1:replicas, inner=(length(D)*n*n)),
                        D=repeat(D, outer=replicas, inner=(n*n)),
                        patch_no=repeat(1:(n*n), outer=(replicas*length(D))),
                        state=Vector{Union{Missing, Int32}}(missing, replicas*length(D)*n*n))

for r=1:replicas

# Import simulation results
df_patch_state = CSV.read(string("../../Output_",version,"/",networkName,"/",networkType,"_out_",networkName,"_patch_state_r",r,".csv"), DataFrame)

# Store in dataframe
patch_state[(patch_state.rep_no.===r),"state"] .= df_patch_state[:,"state"]

for k=1:(length(D)-1)

# Import simulation results
df_p_plants = CSV.read(string("../../Output_",version,"/",networkName,"/",networkType,"_out_",networkName,"_p_plants_r",r,"_D",floor(Int,D[k]*100),".csv"), DataFrame)
df_p_animals = CSV.read(string("../../Output_",version,"/",networkName,"/",networkType,"_out_",networkName,"_p_animals_r",r,"_D",floor(Int,D[k]*100),".csv"), DataFrame)
df_z_plants = CSV.read(string("../../Output_",version,"/",networkName,"/",networkType,"_out_",networkName,"_z_plants_r",r,"_D",floor(Int,D[k]*100),".csv"), DataFrame)
df_z_animals = CSV.read(string("../../Output_",version,"/",networkName,"/",networkType,"_out_",networkName,"_z_animals_r",r,"_D",floor(Int,D[k]*100),".csv"), DataFrame)
df_ab_traits_plants = CSV.read(string("../../Output_",version,"/",networkName,"/",networkType,"_out_",networkName,"_ab_traits_plants_r",r,"_D",floor(Int,D[k]*100),".csv"), DataFrame)
df_ab_traits_animals = CSV.read(string("../../Output_",version,"/",networkName,"/",networkType,"_out_",networkName,"_ab_traits_animals_r",r,"_D",floor(Int,D[k]*100),".csv"), DataFrame)

# Store in dataframe
species_presence[(species_presence.rep_no.===r) .& (species_presence.D.===D[k]) .& (species_presence.level.==="resources"),"p"] .= df_p_plants[:,"p"]
species_presence[(species_presence.rep_no.===r) .& (species_presence.D.===D[k]) .& (species_presence.level.==="consumers"),"p"] .= df_p_animals[:,"p"]

z_values = vcat(z_values, DataFrame(rep_no=fill(r, (size(df_z_plants, 1)+size(df_z_animals, 1))),
                                    D=fill(D[k], (size(df_z_plants, 1)+size(df_z_animals, 1))),
                                    patch_no=vcat(df_z_plants[:,"patch_no"], df_z_animals[:,"patch_no"]),
                                    species=vcat(df_z_plants[:,"species"], df_z_animals[:,"species"]),
                                    level=vcat(fill("resources", size(df_z_plants, 1)), fill("consumers", size(df_z_animals, 1))),
                                    z=vcat(df_z_plants[:,"z"], df_z_animals[:,"z"])))

ab_traits_dt[(ab_traits_dt.rep_no.===r) .& (ab_traits_dt.D.===D[k]) .& (ab_traits_dt.level.==="resources"), "ab"] .= df_ab_traits_plants[:,"ab"]
ab_traits_dt[(ab_traits_dt.rep_no.===r) .& (ab_traits_dt.D.===D[k]) .& (ab_traits_dt.level.==="resources"), "z_mean"] .= df_ab_traits_plants[:,"z_mean"]
ab_traits_dt[(ab_traits_dt.rep_no.===r) .& (ab_traits_dt.D.===D[k]) .& (ab_traits_dt.level.==="resources"), "z_sd"] .= df_ab_traits_plants[:,"z_sd"]
ab_traits_dt[(ab_traits_dt.rep_no.===r) .& (ab_traits_dt.D.===D[k]) .& (ab_traits_dt.level.==="consumers"), "ab"] .= df_ab_traits_animals[:,"ab"]
ab_traits_dt[(ab_traits_dt.rep_no.===r) .& (ab_traits_dt.D.===D[k]) .& (ab_traits_dt.level.==="consumers"), "z_mean"] .= df_ab_traits_animals[:,"z_mean"]
ab_traits_dt[(ab_traits_dt.rep_no.===r) .& (ab_traits_dt.D.===D[k]) .& (ab_traits_dt.level.==="consumers"), "z_sd"] .= df_ab_traits_animals[:,"z_sd"]

end

end

# Write out results
CSV.write(string("../../Results_",version,"/",networkName,"/p_eq.csv"), species_presence)
CSV.write(string("../../Results_",version,"/",networkName,"/z_eq.csv"), z_values)
CSV.write(string("../../Results_",version,"/",networkName,"/ab_traits_dt.csv"), ab_traits_dt)
CSV.write(string("../../Results_",version,"/",networkName,"/patch_state.csv"), patch_state)

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
pmap(run_parallel, ["M_SD_012"], repeat(["mutualistic"], 1), repeat(["v7"], 1), repeat([5], 1))
#pmap(run_parallel, ["A_PH_004"], repeat(["antagonistic"], 1), repeat(["v7"], 1), repeat([5], 1))

#################
# TO RUN THIS SCRIPT
# REPLACE NETWORK NAMES, TYPES AND VERSION IN LINE 126
# DO THE FOLLOWING IN THE COMMAND LINE
# julia -p 5 combine_results_replicas.jl
