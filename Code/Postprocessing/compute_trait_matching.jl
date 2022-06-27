# COMPUTE TRAIT MATCHING

# - this script calculates trait matching
# - for models: v7

using DataFrames, StatsBase, Random, Distributions, CSV, Dates, Distributed, DelimitedFiles, LinearAlgebra

networkName = "M_PL_059"
networkType = "mutualistic"
version = "v7"
n = 100
dD = 0.1


function compute_trait_matching(networkName, networkType, version, n, dD)

# network incidence matrix
M_inc = readdlm(string("../../Data/",networkName,"/M_inc.csv"), ' ', Int)

# number of plant and animas species
n_p = size(M_inc, 1)
n_a = size(M_inc, 2)

# habitat loss
D = collect(0:dD:1)

# theta values
coevolution = readdlm(string("../../Data/",networkName,"/coevolution.csv"), ',', Float64)
theta_values = DataFrame(level=vcat(fill("plants", n_p), fill("animals", n_a)),
                         species=vcat(collect(1:n_p), collect(1:n_a)),
                         theta=coevolution[:,3])

# Create empty dataframe
matching_network = DataFrame(network = fill(networkName, 5*length(D)*n*n),
                             rep_no = repeat(1:5, inner=length(D)*n*n),
                             D=repeat(D, outer=5, inner=(n*n)),
                             patch_no=repeat(1:(n*n), outer=5*length(D)),
                             match=Vector{Union{Missing, Float64}}(missing, 5*length(D)*n*n),
                             match_theta=Vector{Union{Missing, Float64}}(missing, 5*length(D)*n*n),
                             match_resources=Vector{Union{Missing, Float64}}(missing, 5*length(D)*n*n),
                             match_consumers=Vector{Union{Missing, Float64}}(missing, 5*length(D)*n*n))

for r=1:5

for k=1:(length(D)-1)

# Import simulation results
result_plants = CSV.read(string("../../Output_",version,"/",networkName,"/",networkType,"_out_",networkName,"_z_plants_r",r,"_D",floor(Int,D[k]*100),".csv"), DataFrame)
result_animals = CSV.read(string("../../Output_",version,"/",networkName,"/",networkType,"_out_",networkName,"_z_animals_r",r,"_D",floor(Int,D[k]*100),".csv"), DataFrame)

# Patches with data
patches = unique(vcat(result_plants.patch_no, result_animals.patch_no))

for patch=patches

# Filter data for current patch
current_plants = result_plants[result_plants.patch_no .=== patch,:]
current_plants = sort!(current_plants, :species)
current_animals = result_animals[result_animals.patch_no .=== patch,:]
current_animals = sort!(current_animals, :species)

# Plants and animals present
plant_sp = current_plants.species
animal_sp = current_animals.species
n_plants = length(plant_sp)
n_animals = length(animal_sp)

# Define the final trait value of the species after coevolution
z = vcat(current_plants.z, current_animals.z)

theta = vcat(theta_values[(theta_values.level.==="plants") .& in.(theta_values.species, (plant_sp,)),"theta"],
             theta_values[(theta_values.level.==="animals") .& in.(theta_values.species, (animal_sp,)),"theta"])

# Define the sensitivity parameter
alpha = 0.2

# Current incidence matrix
current_incidence = M_inc[plant_sp, animal_sp]

# Calculate trait matching if at least one interaction present
if sum(current_incidence) != 0

# Convert to adjacency matrix
current_adj = hcat(vcat(zeros(Int, n_plants, n_plants), transpose(current_incidence)),
                   vcat(current_incidence, zeros(Int, n_animals, n_animals)))

# Assign NA to
current_adj = convert(Array{Union{Missing,Int64}}, current_adj)
current_adj[current_adj.===0] .= missing

# Calculate the degree of trait matching between species of different guilds
matching_diff = [norm(z[i]-z[j]) for i in eachindex(z), j in eachindex(z)]
matching = exp.(-alpha.*(matching_diff).^2)

# Guild trait matching
plants_matching = (sum(matching[1:n_plants, 1:n_plants]) - n_plants) / (n_plants*n_plants-n_plants)
animals_matching = (sum(matching[(n_plants+1):end, (n_plants+1):end]) - n_animals) / (n_animals*n_animals-n_animals)

# Assign NA to matching between non-interacting species
matching = matching .* current_adj

# Calculate the mean trait matching value for the network
network_matching = mean(skipmissing(matching))

# Store results in dataframe
matching_network[(matching_network.rep_no.===r) .& (matching_network.D.===D[k]) .& (matching_network.patch_no.===patch), "match"] .= network_matching
matching_network[(matching_network.rep_no.===r) .& (matching_network.D.===D[k]) .& (matching_network.patch_no.===patch), "match_resources"] .= plants_matching
matching_network[(matching_network.rep_no.===r) .& (matching_network.D.===D[k]) .& (matching_network.patch_no.===patch), "match_consumers"] .= animals_matching

end

# Calculate the degree of trait matching with theta
matching_theta = exp.(-alpha*(theta - z).^2)

# Calculate the mean trait matching value for the network
network_matching_theta = mean(skipmissing(matching_theta))

# Store results in dataframe
matching_network[(matching_network.rep_no.===r) .& (matching_network.D.===D[k]) .& (matching_network.patch_no.===patch), "match_theta"] .= network_matching_theta

end
end
end

# Writeout results
CSV.write(string("../../Results_",version,"/",networkName,"/matching_network.csv"), matching_network)

return 0
end


@time compute_trait_matching(networkName, networkType, version, n, dD)
