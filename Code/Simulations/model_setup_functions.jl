# functions used to setup simulations

function setup_habitat_destruction(dD, n)

	# number of patches
	n_patch = n * n

	# number of patches destroyed in a step
	dn = floor(Int, dD*n_patch)

	# vector of number of destroyed patches
	d = collect(0:dn:n_patch)

	# vector of fraction of destroyed patches
	D = d/n_patch

	# return
	return n_patch, dn, d, D

end

function setup_grids(n, n_p, n_a, d, plant_parameters, animal_parameters)

	# grid of patch states (1=pristine, 0=destroyed) for each D
	x_state = ones(Int, n, n, length(d))

	# grid of plants (1=present, 0=absent) for each species
	x_p = ones(Int, n, n, n_p)

	# grid of animals (1=present, 0=absent) for each species
	x_a = ones(Int, n, n, n_a)

	# grid of plant trait values
	z_p = Array{Union{Missing, Float64}}(missing, n, n, n_p)

	# grid of animal trait values
	z_a = Array{Union{Missing, Float64}}(missing, n, n, n_a)

	# set initial trait values
	for h=1:n_p
		# initially trait value = theta
    	z_p[:,:,h] .= plant_parameters[3,h]
  	end

  	for h=1:n_a
		# initially trait value = theta
    	z_a[:,:,h] .= animal_parameters[3,h]
  	end

	# return
	return x_state, x_p, x_a, z_p, z_a

end


function initialise_dataframes_store_results(tmax, n)

	# initialise dataframes for storing results
	df_ab_traits = DataFrame(t = repeat(1:tmax, inner=n),
	                         species = repeat(1:n, outer=tmax),
	                         ab = Vector{Union{Missing, Float32}}(missing, tmax*n),
	                         z_mean = Vector{Union{Missing, Float32}}(missing, tmax*n),
	                         z_sd = Vector{Union{Missing, Float32}}(missing, tmax*n))

	# return
	return df_ab_traits

end

function setup_moore_neighborhood(i, j, n)

	# Moore's neighborhood
	neigh = [i-1 j-1
             i+0 j-1
             i+1 j-1
             i-1 j+0
	         i+1 j+0
    	     i-1 j+1
             i+0 j+1
             i+1 j+1]

	# remove non-existent neighbours
	neigh = neigh[(neigh[:,1].>=1) .& (neigh[:,2].>=1) .& (neigh[:,1].<=n) .&(neigh[:,2].<=n),:]

	# return
	return neigh

end

function setup_coevolution(current_inc, p_current_sp, a_current_sp, p_nonint_sp, a_nonint_sp, M_inc, plant_parameters, animal_parameters, z_p, z_a, i, j)

	# remove empty rows and columns (non-interacting speices)
    p_int = vec(mapslices(col -> any(col .!== 0), current_inc, dims = 2))
	a_int = vec(mapslices(row -> any(row .!== 0), current_inc, dims = 1))
	p_int_sp = p_current_sp[p_int]
    a_int_sp = a_current_sp[a_int]

    # non-interacting species
    p_nonint_sp = vcat(p_nonint_sp, setdiff(p_current_sp, p_int_sp))
    a_nonint_sp = vcat(a_nonint_sp, setdiff(a_current_sp, a_int_sp))

    # current incidence matrix
    current_inc_coev = M_inc[p_int_sp, a_int_sp]

    # current number of plants and animals
    n_p_current = length(p_int_sp)
    n_a_current = length(a_int_sp)
    n_sp_current = n_p_current + n_a_current

    # build adjacency matrix
    f = vcat(hcat(zeros(Int, size(current_inc_coev,1),size(current_inc_coev,1)), current_inc_coev),hcat(transpose(current_inc_coev), zeros(Int,size(current_inc_coev,2),size(current_inc_coev,2))))

    # coevolution parameters
    m = vcat(plant_parameters[1, p_int_sp], animal_parameters[1, a_int_sp])
    phi = vcat(plant_parameters[2, p_int_sp], animal_parameters[2, a_int_sp])
    theta = vcat(plant_parameters[3, p_int_sp], animal_parameters[3, a_int_sp])

    # trait values
    z = vcat(z_p[i,j,p_int_sp], z_a[i,j,a_int_sp])

	# return
	return p_int_sp, a_int_sp, p_nonint_sp, a_nonint_sp, n_p_current, n_a_current, f, m, phi, theta, z

end
