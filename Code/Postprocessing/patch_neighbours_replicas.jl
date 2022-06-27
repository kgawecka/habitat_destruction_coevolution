# PATCH NEIGHBOURS

# - this script determines the number of pristine neighbours of each patch
# - for models: v7, v7a

@everywhere using DataFrames, StatsBase, Random, Distributions, CSV, Dates, Distributed, DelimitedFiles


@everywhere function patch_neighbours(networkName, networkType, version, n, dD)

  # habitat loss
  D = collect(0:dD:1)

  patch_neigh = DataFrame(network = fill(networkName, 5*length(D)*n*n),
                          rep_no = repeat(1:5, inner=length(D)*n*n),
                          D = repeat(D, outer=5, inner=n*n),
                          patch_no = repeat(1:(n*n), outer=5*length(D)),
                          no_neigh = Vector{Union{Missing, Int64}}(missing, 5*length(D)*n*n),     # number of pristine neighbours
                          no_occ_neigh = Vector{Union{Missing, Int64}}(missing, 5*length(D)*n*n), # number of occupied neighbours
                          no_species = Vector{Union{Missing, Int64}}(missing, 5*length(D)*n*n))   # total number of neighbouring species

  # import presnece/absence data
  p_eq = CSV.read(string("../../Results_",version,"/",networkName,"/p_eq.csv"), DataFrame)
  p_eq = dropmissing(p_eq)

  # import patch state data
  state_df = CSV.read(string("../../Results_",version,"/",networkName,"/patch_networks.csv"), DataFrame)
  state_df = state_df[:,["rep_no", "D", "patch_no", "state"]]

  # import patch coordinates
  patchXY = CSV.read(string("../../Data/patchXY.csv"), DataFrame)

  for r=1:5

    for k=1:length(D)

      # data for current D
      data = p_eq[(p_eq.rep_no.===r) .& (p_eq.D.===D[k]),:]

      # occupied patches
      patches = unique(data[data.p.===1,"patch_no"])

      if length(patches)===0
        break
      end

      for i=1:length(patches)

        current_patch = patches[i]
        currentX = patchXY[current_patch,"X"]
        currentY = patchXY[current_patch,"Y"]

        # Moore's neighborhood
        neigh = [currentX-1 currentY-1
                 currentX+0 currentY-1
                 currentX+1 currentY-1
                 currentX-1 currentY+0
                 currentX+1 currentY+0
                 currentX-1 currentY+1
                 currentX+0 currentY+1
                 currentX+1 currentY+1]
        neigh = neigh[(neigh[:,1].>=1) .& (neigh[:,2].>=1) .& (neigh[:,1].<=n) .& (neigh[:,2].<=n),:]   # remove non-existent neighbours

        neigh = DataFrame(X=neigh[:,1], Y=neigh[:,2])
        neigh = leftjoin(neigh, patchXY, on=["X", "Y"])

        neigh_no = neigh.patch_no
        neigh_occ = data[in.(data.patch_no, (neigh_no,)) .& (data.p.===1),:]

        neigh_pristine = state_df[(state_df.rep_no.===r) .& (state_df.D.===D[k]) .& in.(state_df.patch_no, (neigh_no,)) .& (state_df.state.===1),:]

        # store results
        patch_neigh[(patch_neigh.rep_no.===r) .& (patch_neigh.D.===D[k]) .& (patch_neigh.patch_no.===current_patch),"no_neigh"] .= length(neigh_pristine.patch_no)
        patch_neigh[(patch_neigh.rep_no.===r) .& (patch_neigh.D.===D[k]) .& (patch_neigh.patch_no.===current_patch),"no_occ_neigh"] .= length(unique(neigh_occ.patch_no))
        patch_neigh[(patch_neigh.rep_no.===r) .& (patch_neigh.D.===D[k]) .& (patch_neigh.patch_no.===current_patch),"no_species"] .= length(neigh_occ.species)

      end # i loop
    end # k loop
  end # r loop

  # write out results
  CSV.write(string("../../Results_",version,"/",networkName,"/patch_neighbours.csv"), patch_neigh)

  return 0
end


# Function for running networks in parallel
@everywhere function run_parallel(network, type, model_version)

  networkName = network
  networkType = type
  version = model_version
  n = 100
  dD = 0.1

  patch_neighbours(networkName, networkType, version, n, dD)

end


# Run simulations
#pmap(run_parallel, [networks], repeat([type], length(networks)), repeat([version], length(networks)))
pmap(run_parallel, ["M_PL_006", "M_PL_010", "M_PL_036", "M_PL_037", "M_PL_059"], repeat(["mutualistic"], 5), repeat(["v7"], 5))

#################
# TO RUN THIS SCRIPT
# REPLACE NETWORK NAMES, TYPES AND VERSION IN LINE 183
# DO THE FOLLOWING IN THE COMMAND LINE
# julia -p 5 patch_neighbours_replicas.jl
