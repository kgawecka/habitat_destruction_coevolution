# GENERATE SPECIES PARAMETERS - MUTUALISM

using DataFrames, StatsBase, Random, Distributions, CSV, Dates, Distributed, DelimitedFiles, LinearAlgebra

networkName = "M_PL_006"

# set working directory
cd(string("../../Data/",networkName,"/"))

# network incidence matrix
M_inc = readdlm("M_inc.csv", ' ', Int)

# number of plants and animals
n_p = size(M_inc, 1)
n_a = size(M_inc, 2)
n_sp = n_p + n_a


# extinction and colonisation probabilities

ep = rand(Uniform(0.14, 0.16), n_p)
ea = rand(Uniform(0.14, 0.16), n_a)
cp = rand(Uniform(0.04, 0.06), n_p, n_a) .* M_inc
ca = rand(Uniform(0.04, 0.06), n_p, n_a) .* M_inc

writedlm(string("ep.csv"), ep, ",")
writedlm(string("ea.csv"), ea, ",")
writedlm(string("cp.csv"), cp, ",")
writedlm(string("ca.csv"), ca, ",")


# coevolution parameters

# Sample phi values
phi = rand(Normal(0.5, 0.01), n_sp)
# Make sure no values are negative
while any(phi.<0)
phi = rand(Normal(0.5, 0.01), n_sp)
end

# Sample m values
m = rand(Normal(0.7, 0.01), n_sp)
# Make sure no values are negative
while any(m.<0)
m = rand(Normal(0.7, 0.01), n_sp)
end

coevolution = hcat(m,
                   phi,
                   rand(Uniform(0, 10), n_sp))  # theta

writedlm(string("coevolution.csv"), coevolution, ",")
