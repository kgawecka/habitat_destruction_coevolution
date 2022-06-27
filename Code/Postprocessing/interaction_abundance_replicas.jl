# INTERACTION ABUNDANCES

# - this script calculates interaction abundances
# - for models: v7, v7a

# Packages
@everywhere using DataFrames, StatsBase, Random, Distributions, CSV, Dates, DelimitedFiles, LinearAlgebra


@everywhere function interaction_abundance(networkName, version, replicas)

n = 100

# network incidence matrix
M_inc = readdlm(string("../../Data/",networkName,"/M_inc.csv"), ',', Int)

# number of plant and animal species, and interactions
n_p = size(M_inc, 1)
n_a = size(M_inc, 2)
n_int = sum(M_inc)

# import presence/absence data
p_eq = CSV.read(string("../../Results_",version,"/",networkName,"/p_eq.csv"), DataFrame)

# habitat loss
D = unique(p_eq.D)

# create empty dataframe
abundance_int = DataFrame(rep_no = repeat(1:replicas, inner=length(D)*n_int),
                          D = repeat(D, outer=replicas, inner=n_int),
                          int_no = repeat(1:n_int, outer=replicas*length(D)),
                          abundance = fill(0.0, replicas*length(D)*n_int))

for r=1:replicas
for k=1:length(D)

# data for current D
p_data = p_eq[(p_eq.rep_no.===r) .& (p_eq.D.===D[k]),:]

# interaction number count
int = 0

for i=1:n_p
for j=1:n_a

if M_inc[i,j] === 0
continue
else
  # update interaction number count
  int = int + 1

# patches where plant i in present
patches_plant = p_data[(p_data.level.==="resources") .& (p_data.species.===i) .& (p_data.p.===1),"patch_no"]
# patches where animal j in present
patches_animal = p_data[(p_data.level.==="consumers") .& (p_data.species.===j) .& (p_data.p.===1),"patch_no"]

# patches with both species present
patches_both = intersect(patches_plant, patches_animal)

# store abundance
abundance_int[(abundance_int.rep_no.===r) .& (abundance_int.D.===D[k]) .& (abundance_int.int_no.===int),"abundance"] .= length(patches_both)/(n*n)
end

end #j loop
end # i loop
end # k loop
end # r loop

# writeout results
CSV.write(string("../../Results_",version,"/",networkName,"/abundance_interactions.csv"), abundance_int)

return 0
end


# Function for running networks in parallel
@everywhere function run_parallel(network, model_version, reps)

networkName = network
version = model_version
replicas = reps

interaction_abundance(networkName, version, replicas)

end


# Run simulations
#pmap(run_parallel, [networks], repeat([version], length(networks)), repeat([replicas], length(networks)))
pmap(run_parallel, ["M_PL_006", "M_PL_010", "M_PL_036", "M_PL_037", "M_PL_059"], repeat(["v7"], 5), repeat([1], 5))


#################
# TO RUN THIS SCRIPT
# REPLACE NETWORK NAMES, TYPES AND VERSION IN LINE 184
# DO THE FOLLOWING IN THE COMMAND LINE
# julia -p 5 interaction_abundance_replicas.jl
