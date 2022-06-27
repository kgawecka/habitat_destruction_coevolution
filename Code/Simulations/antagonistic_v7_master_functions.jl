function simulate_coevolution_or_evolution(current_inc, p_current_sp, a_current_sp, p_nonint_sp, a_nonint_sp, M_inc, plant_parameters, animal_parameters, z_p, i, j, alpha, sigma, z_a)

	# check if at least one interaction present -> coevolution
    	if sum(current_inc) !== 0

    	# setup parameters for coevolution
		p_int_sp, a_int_sp, p_nonint_sp, a_nonint_sp, n_p_current, n_a_current, f, m, phi, theta, z = setup_coevolution(current_inc, p_current_sp, a_current_sp, p_nonint_sp, a_nonint_sp, M_inc, plant_parameters, animal_parameters, z_p, z_a, i, j)

		# simulate coevolution
		z_p, z_a, q_m = simulate_coevolution(f, z, alpha, m, n_p_current, n_a_current, sigma, phi, z_p, z_a, theta, i, j, p_int_sp, a_int_sp)

    end

    # if non-interacting species present -> environmental selection only
    if length(p_nonint_sp) !== 0 || length(a_nonint_sp) !== 0

    # simulate evolution (only environmental selection)
    z_p, z_a = simulate_evolution(plant_parameters, animal_parameters, p_nonint_sp, a_nonint_sp, z_p, z_a, i, j)

    end

	# return
	return z_p, z_a

end


function update_patch_sequentially(n, k, x_state, z_p, z_a, n_p, x_p_old, M_inc, x_a_old, z_p_old, z_a_old, ep0, ep, x_p, cp, n_a, ea, x_a, alpha, ca, plant_parameters, animal_parameters, sigma)

	# update each patch sequentially
	for j=1:n
    	for i=1:n

    	# check state of the current patch
        if x_state[i,j,k] === 0
        	z_p[i,j,:] .= missing
            z_a[i,j,:] .= missing
            # if patch is destroyed, move to next patch
            continue
        end

        # setup moore neighborhood
        neigh = setup_moore_neighborhood(i, j, n)

    	# resource extinctions and colonisations
        x_p, z_p = resource_extinctions_and_colonisations(n_p, x_p_old, i, j, M_inc, x_a_old, z_p_old, z_a_old, ep0, ep, x_p, z_p, neigh, cp, alpha, sigma)

        # consumer extinctions and colonisations
        x_a, z_a = consumer_extinctions_and_colonisations(n_a, x_a_old, i, j, ea, x_a, z_a, M_inc, x_p_old, z_p_old, neigh, z_a_old, alpha, ca)

        # species present in current patch
    	p_current_sp = findall(x_p[i,j,:].===1)
        a_current_sp = findall(x_a[i,j,:].===1)

        # check if at least one species present
        if length(p_current_sp) === 0 && length(a_current_sp) === 0
        	continue
        end

        # update grid after colonisation extinction
        current_inc, p_nonint_sp, a_nonint_sp = update_grid_after_colonisation_extinction(M_inc, p_current_sp, a_current_sp)

        # simulate coevolution or evolution
        z_p, z_a = simulate_coevolution_or_evolution(current_inc, p_current_sp, a_current_sp, p_nonint_sp, a_nonint_sp, M_inc, plant_parameters, animal_parameters, z_p, i, j, alpha, sigma, z_a)

        end # i loop
    end # j loop

	# return
	return x_p, z_p, x_a, z_a
end

function store_and_write_results(n_p, x_p, n_a, x_a, z_p, z_a, t_sim, n_patch, df_ab_traits_p, df_ab_traits_a, networkName, rep_no, D, k)

	# store and output results
	df_p_plants = build_p_data(n_p, x_p, t_sim, n_patch)
	df_p_animals = build_p_data(n_a, x_a, t_sim, n_patch)
	df_z_plants = build_z_data(n_p, z_p, t_sim, n_patch)
	df_z_animals = build_z_data(n_a, z_a, t_sim, n_patch)

	# output results
    CSV.write(string("../../Output_v7/",networkName,"/antagonistic_out_",networkName,"_p_plants_r",rep_no,"_D",floor(Int,D[k]*100),".csv"), df_p_plants)
    CSV.write(string("../../Output_v7/",networkName,"/antagonistic_out_",networkName,"_p_animals_r",rep_no,"_D",floor(Int,D[k]*100),".csv"), df_p_animals)
    CSV.write(string("../../Output_v7/",networkName,"/antagonistic_out_",networkName,"_z_plants_r",rep_no,"_D",floor(Int,D[k]*100),".csv"), df_z_plants)
    CSV.write(string("../../Output_v7/",networkName,"/antagonistic_out_",networkName,"_z_animals_r",rep_no,"_D",floor(Int,D[k]*100),".csv"), df_z_animals)
    CSV.write(string("../../Output_v7/",networkName,"/antagonistic_out_",networkName,"_ab_traits_plants_r",rep_no,"_D",floor(Int,D[k]*100),".csv"), df_ab_traits_p)
	CSV.write(string("../../Output_v7/",networkName,"/antagonistic_out_",networkName,"_ab_traits_animals_r",rep_no,"_D",floor(Int,D[k]*100),".csv"), df_ab_traits_a)

	# return
	return 1

end



function iterate_model_until_steady_state(tmax, x_p, x_a, z_p, z_a, n, n_p, n_a, k, x_state, M_inc, ep0, ep, cp, ea, alpha, ca, plant_parameters, animal_parameters, sigma, df_ab_traits_p, df_ab_traits_a, n_patch, tol, tmin, t_sim)

	# iterate until state
    for g=2:(tmax+1)

    	# make copies of variables
        x_p_old = copy(x_p)
        x_a_old = copy(x_a)
        z_p_old = copy(z_p)
        z_a_old = copy(z_a)

        # update each patch sequentially
        x_p, z_p, x_a, z_a = update_patch_sequentially(n, k, x_state, z_p, z_a, n_p, x_p_old, M_inc, x_a_old, z_p_old, z_a_old, ep0, ep, x_p, cp, n_a, ea, x_a, alpha, ca, plant_parameters, animal_parameters, sigma)

        # calculate global abundance and store results
        df_ab_traits_p = calculate_global_abundance_and_traits_and_store(x_p, z_p, n_patch, n_p, df_ab_traits_p, t_sim)
        df_ab_traits_a = calculate_global_abundance_and_traits_and_store(x_a, z_a, n_patch, n_a, df_ab_traits_a, t_sim)

        t_sim = g

        # check convergence
        if g > tmin

            # check model convergence
            dp, da = check_model_convergence(x_p_old, n_patch, n_p, x_p, x_a_old, n_a, x_a)
            if all(i -> i < tol, dp) && all(i -> i < tol, da)
                break
            end
        end

      end # g loop

	# return
	return n_p, x_p, n_a, x_a, z_p, z_a, n_patch, df_ab_traits_p, df_ab_traits_a, k


end
