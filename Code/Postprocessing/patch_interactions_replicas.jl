# PATCH INTERACTIONS, SPECIES, CONNECTANCE, NESTEDNESS

# - this script calculates size, connectance and nestedness of local networks
# - for models: v7, v7a

@everywhere using DataFrames, StatsBase, Random, Distributions, CSV, Dates, Distributed, DelimitedFiles


@everywhere function nestedness_Fortuna(B)

    # nestedness of rows
    nested_rows = 0
    for i=1:size(B,1)
        for j=1:size(B,1)
            if j>i
                shared=sum(B[i,:].*B[j,:]) # sum of common interactions
                min_shared = min(sum(B[i,:]),sum(B[j,:])) # min of the degrees
                nested_rows = nested_rows+(shared/min_shared)
            end
        end
    end

    # nestedness of columns
    nested_columns = 0
    for i=1:size(B,2)
        for j=1:size(B,2)
            if j>i
                shared=sum(B[:,i].*B[:,j]) # sum of common interactions
                min_shared = min(sum(B[:,i]),sum(B[:,j])) # min of the degrees
                nested_columns = nested_columns+(shared/min_shared)
            end
        end
    end

    # nestedness of the network
    nestedness_network = (nested_rows+nested_columns)/((size(B,1)*(size(B,1)-1)/2)+(size(B,2)*(size(B,2)-1)/2))

    return nestedness_network
end


@everywhere function patch_networks(networkName, networkType, version, n, dD)

# network incidence matrix
M_inc = readdlm(string("../../Data/",networkName,"/M_inc.csv"), ' ', Int)

# number of plants, animals and interactions
n_p = size(M_inc, 1)
n_a = size(M_inc, 2)
n_int = sum(M_inc)

# habitat loss vector
D = collect(0:dD:1)

# number of patches
n_patch = n*n

# initialise dataframe for storing results
df_int = DataFrame(network = fill(networkName, 5*n_patch*length(D)),
                   rep_no = repeat(1:5, inner=n_patch*length(D)),
                   D = repeat(D, outer=5, inner=n_patch),
                   patch_no = repeat(1:n_patch, outer=5*length(D)),
                   state = Vector{Union{Missing, Int64}}(missing, 5*n_patch*length(D)),  # patch state
                   int = Vector{Union{Missing, Float64}}(missing, 5*n_patch*length(D)),  # fraction of interactions
                   con = Vector{Union{Missing, Float64}}(missing, 5*n_patch*length(D)),  # connectance
                   nest = Vector{Union{Missing, Float64}}(missing, 5*n_patch*length(D)), # nestedness
                   n_p = Vector{Union{Missing, Int64}}(missing, 5*n_patch*length(D)),    # no of plants
                   n_a = Vector{Union{Missing, Int64}}(missing, 5*n_patch*length(D)))    # no of animals

for r=1:5

# import patch state data
df_state = CSV.read(string("../../Output_",version,"/",networkName,"/",networkType,"_out_",networkName,"_patch_state_r",r,".csv"), DataFrame)

df_int[(df_int.rep_no.===r),"state"] .= df_state.state

for k=1:(length(D)-1)

# import simulation results
df_p_plants = CSV.read(string("../../Output_",version,"/",networkName,"/",networkType,"_out_",networkName,"_p_plants_r",r,"_D",floor(Int,D[k]*100),".csv"), DataFrame)
df_p_animals = CSV.read(string("../../Output_",version,"/",networkName,"/",networkType,"_out_",networkName,"_p_animals_r",r,"_D",floor(Int,D[k]*100),".csv"), DataFrame)

for p=1:n_patch

# if patch is destroyed, move to next
if df_state[(df_state.D.===D[k]) .& (df_state.patch_no.===p), "state"] === 0
continue
end

# data for current patch
p_current = df_p_plants[df_p_plants.patch_no.===p,"p"]
a_current = df_p_animals[df_p_animals.patch_no.===p,"p"]
n_p_current = sum(p_current)
n_a_current = sum(a_current)

# if patch is empty, move to next
if n_p_current === 0 && n_a_current === 0
continue
end

# current incidence matrix
M_patch = ones(Int, (n_p,n_a)) .* p_current .* transpose(a_current) .* M_inc

# species
df_int[(df_int.rep_no.===r) .& (df_int.D.===D[k]) .& (df_int.patch_no.===p),"n_p"] .= n_p_current
df_int[(df_int.rep_no.===r) .& (df_int.D.===D[k]) .& (df_int.patch_no.===p),"n_a"] .= n_a_current

# interactions
int_patch = sum(M_patch)
df_int[(df_int.rep_no.===r) .& (df_int.D.===D[k]) .& (df_int.patch_no.===p),"int"] .= int_patch / n_int

# if at least one plant and one animal present
if n_p_current > 0 && n_a_current > 0
# connectance
con_patch = int_patch / (n_p_current * n_a_current)
df_int[(df_int.rep_no.===r) .& (df_int.D.===D[k]) .& (df_int.patch_no.===p),"con"] .= con_patch

# remove empty rows or columns
M_patch = M_patch[vec(sum(M_patch, dims=2).>0),vec(sum(M_patch, dims=1).>0)]

# if at least two plants and two animals present
if size(M_patch,1) > 1 && size(M_patch,2) > 1

# sort rows and columns by decreasing number of links
rowsums = DataFrame(rows=collect(1:size(M_patch,1)), sums=sum(M_patch, dims=2)[:,1])
colsums = DataFrame(cols=collect(1:size(M_patch,2)), sums=sum(M_patch, dims=1)[1,:])
sort!(rowsums, :sums, rev=true)
sort!(colsums, :sums, rev=true)
M_patch = M_patch[rowsums.rows,colsums.cols]

# nestedness using Fortuna et al. (2019) definition
df_int[(df_int.rep_no.===r) .& (df_int.D.===D[k]) .& (df_int.patch_no.===p),"nest"] .= nestedness_Fortuna(M_patch)

end
end

end # p loop
end # k loop
end # r loop

# remove empty rows
#df_int = dropmissing(df_int)

# write out results
CSV.write(string("../../Results_",version,"/",networkName,"/patch_networks.csv"), df_int)

return 0
end


# Function for running networks in parallel
@everywhere function run_parallel(network, type, model_version)

  networkName = network
  networkType = type
  version = model_version
  n = 100
  dD = 0.1

  patch_networks(networkName, networkType, version, n, dD)

end


# Run simulations
#pmap(run_parallel, [networks], repeat([type], length(networks)), repeat([version], length(networks)))
pmap(run_parallel, ["M_PL_006", "M_PL_010", "M_PL_036", "M_PL_037", "M_PL_059"], repeat(["mutualistic"], 5), repeat(["v7"], 5))

#################
# TO RUN THIS SCRIPT
# REPLACE NETWORK NAMES, TYPES AND VERSION IN LINE 183
# DO THE FOLLOWING IN THE COMMAND LINE
# julia -p 5 patch_interactions_replicas.jl
