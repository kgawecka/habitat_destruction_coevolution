# load packages
using Distributed, DataFrames, StatsBase, Random, Distributions, CSV, Dates, DelimitedFiles, LinearAlgebra
@everywhere using DataFrames, StatsBase, Random, Distributions, CSV, Dates, DelimitedFiles, LinearAlgebra

# load functions
@everywhere include("model_setup_functions.jl")
@everywhere include("coevolution_auxiliary_functions.jl")
@everywhere include("construct_output_files_functions.jl")
@everywhere include("random_habitat_destruction_functions.jl")
@everywhere include("mutualistic_extinction_colonisation_functions.jl")
@everywhere include("mutualistic_v7a_master_functions.jl")

# Define model function
@everywhere function dynamicsD(n_p, n_a, n, dD, tmin, tmax, tol, ep, ea, cp, ca, M_inc, rep_no, networkName, plant_parameters, animal_parameters, alpha)

    # setup variables for habitat destruction
    n_patch, dn, d, D = setup_habitat_destruction(dD, n)

    # setup grids
    x_state, x_p, x_a, z_p, z_a = setup_grids(n, n_p, n_a, d, plant_parameters, animal_parameters)

    @inbounds for k=1:(length(d)-1)

        # counter
        t_sim = 1

        # run colonisation/extinction and coevolution until steady stated
        n_p, x_p, n_a, x_a, z_p, z_a, n_patch, k = iterate_model_until_steady_state(tmax, x_p, x_a, z_p, z_a, n, n_p, n_a, k, x_state, M_inc, ep, cp, ea, alpha, ca, plant_parameters, animal_parameters, n_patch, tol, tmin)

        # store and output results
        store_results = store_and_write_results(n_p, x_p, n_a, x_a, z_p, z_a, t_sim, n_patch, networkName, rep_no, D, k)

        # habitat destruction
        x_p, x_a, x_state  = habitat_destruction(n_patch, x_state, k, dn, n, x_p, x_a)

    end # k loop

  # store and output results
  df_state = DataFrame(D = repeat(D, inner=n_patch),
                       patch_no = repeat(1:n_patch, outer=length(D)),
                       state = vec(x_state))
  CSV.write(string("../../Output_v7a/",networkName,"/mutualistic_out_",networkName,"_patch_state_r",rep_no,".csv"), df_state)


  # return
  return 0

end


# Define function for running replicas
@everywhere function run_parallel(replicate, networkName)

  # simulation parameters
  n = 100            # grid size nxn
  dD = 0.1           # fraction of patches destoryed at a time
  tmin = 10          # minimum number of timesteps
  tmax = 1000        # maximum number of timesteps
  tol = 0.001        # convergence tolerance

  # network incidence matrix
  M_inc = readdlm(string("../../Data/",networkName,"/M_inc.csv"), ' ', Int)
  # number of plants and animals
  n_p = size(M_inc, 1)
  n_a = size(M_inc, 2)
  n_sp = n_p + n_a

  # extinction and colonisation probabilities
  ep = readdlm(string("../../Data/",networkName,"/ep.csv"), ',', Float64)
  ea = readdlm(string("../../Data/",networkName,"/ea.csv"), ',', Float64)
  cp = readdlm(string("../../Data/",networkName,"/cp.csv"), ',', Float64)
  ca = readdlm(string("../../Data/",networkName,"/ca.csv"), ',', Float64)

  # coevolution parameters
  coevolution = readdlm(string("../../Data/",networkName,"/coevolution.csv"), ',', Float64)
  plant_parameters = transpose(coevolution[1:n_p,:])
  animal_parameters = transpose(coevolution[(n_p+1):end,:])
  alpha = 0.2

  rep_no = replicate

  dynamicsD(n_p, n_a, n, dD, tmin, tmax, tol, ep, ea, cp, ca, M_inc, rep_no, networkName, plant_parameters, animal_parameters, alpha)
end

# Define network to run
networkName = string(ARGS[1])

# Run simulations
pmap(run_parallel, 1:5, repeat([networkName], 5))


#################
# TO RUN THIS SCRIPT
# DO THE FOLLOWING IN THE COMMAND LINE, REPLACE NETWORK NAME
# julia -p 5 fragmentation_mutualistic_coevolution_v7a_modular.jl "M_SD_025"
