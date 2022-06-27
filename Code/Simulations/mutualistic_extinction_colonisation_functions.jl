function get_animal_partners_in_current_patch(i, j, h, M_inc, x_a_old)

	# get animal partners of plant h
	a_part = findall(M_inc[h,:].===1)

	# get animals present in current patch
	a_neigh_all = findall(x_a_old[i,j,:].===1)

	# get animal partners present in current patch
	a_part_sp = a_neigh_all[findall(x-> x in a_part, a_neigh_all)]

	# return
	return a_part_sp

end

function resource_extinctions(h, i, j, x_p, z_p, ep)

	if rand() < ep[h]
		x_p[i,j,h] = 0
		z_p[i,j,h] = missing

	end

	# return grid of plants and grid of plant trait values
	return x_p, z_p

end

function get_animal_partners_in_current_patch(M_inc, h)

	# find animal partners of plant h
	a_part = findall(M_inc[h,:].===1)

	# return
	return a_part

end

function compute_resource_colonisation_prob(z_a_old, z_p_old, i_n, j_n, h, a_part_sp, alpha, cp)

	# trait value of plant
	z_p_s = z_p_old[i_n,j_n,h]

	# trais of animal partners
	z_a_part = z_a_old[i_n,j_n,a_part_sp]

	# trait matching
	p_match = exp.(-alpha.*((z_p_s .- z_a_part).^2))

	# "not colonised" probability
	pc_p = 1
	j_count = 0

	for l=a_part_sp

		j_count = j_count+1
		prod = (1 - cp[h,l] / j_count) * (1 - p_match[j_count] * cp[h,l])
		pc_p = pc_p * prod

	end

	# return
	return pc_p


end

function resource_colonisation_prob(neigh, x_p_old, x_a_old, a_part, z_a_old, z_p_old, h, alpha, cp, i, j, x_p, z_p)

	# check neighbours
    for f=1:size(neigh,1)

		# neighbour coordinates
    	i_n = neigh[f,1]
    	j_n = neigh[f,2]

    	if x_p_old[i_n,j_n,h] === 0
    		continue

    	# if plant h present in neighbouring patch
		else

			# animals present in neighbour patch
			a_neigh_all = findall(x_a_old[i_n,j_n,:].===1)

			# animal partners present in neighbour patch
			a_part_sp = a_neigh_all[findall(x-> x in a_part, a_neigh_all)]

			# not colonised probability
			if length(a_part_sp) != 0 # if plant species h and its partner(s) present in neighbouring patch

				pc_p = compute_resource_colonisation_prob(z_a_old, z_p_old, i_n, j_n, h, a_part_sp, alpha, cp)

				if rand() > pc_p
					x_p[i,j,h] = 1
					z_p[i,j,h] = z_p_old[i_n,j_n,h]
					break   # stop colonisations
				end

			end

    	end

    end

	# return
	return(x_p, z_p)

end

function resource_colonisations(neigh, x_p_old, z_a_old, z_p_old, cp, i, j, h, x_p, z_p, M_inc, x_a_old, alpha)

	# find animal partners of plant
	a_part = get_animal_partners_in_current_patch(M_inc, h)

	x_p, z_p = resource_colonisation_prob(neigh, x_p_old, x_a_old, a_part, z_a_old, z_p_old, h, alpha, cp, i, j, x_p, z_p)

	# return
	return x_p, z_p

end

function consumer_extinctions(i, j, h, ea, x_a, z_a)

	if rand() < ea[h]
		x_a[i,j,h] = 0
		z_a[i,j,h] = missing
	end

	# return
	return x_a, z_a

end

function get_plant_partners_in_current_patch(h, M_inc, x_p_old, i, j)

	# get plant partners of plant h
	p_part = findall(M_inc[:,h].===1)

	# get plants present in current patch
	p_neigh_all = findall(x_p_old[i,j,:].===1)

	# get plant partners present in current patch
	p_part_sp = p_neigh_all[findall(x-> x in p_part, p_neigh_all)]

	# return
	return p_part_sp

end

function compute_consumer_colonisation_prob(z_p_old, i, j, p_part_sp, neigh, x_a_old, h, z_a_old, alpha, ca, x_a, z_a)

	# traits of plant partners
	z_p_part = z_p_old[i,j,p_part_sp]

	for f = 1:size(neigh,1)
		i_n = neigh[f,1]
		j_n = neigh[f,2]

		if x_a_old[i_n, j_n, h] === 0
			continue

		# if animal h present in neighbouring patch
		else
			# trait value of animal
			z_a_s = z_a_old[i_n,j_n,h]
			# trait matching
            p_match = exp.(-alpha.*((z_a_s .- z_p_part).^2))
			# "not colonised" probability
            pc_a = 1
            i_count = 0

            for m=p_part_sp
            	i_count = i_count+1
            	prod = (1 - ca[m,h] / i_count) * (1 - p_match[i_count] * ca[m,h])
    		    pc_a = pc_a * prod
            end

			if rand() > pc_a
            	x_a[i,j,h] = 1
            	z_a[i,j,h] = z_a_old[i_n,j_n,h]
				# stop colonisations
                break
            end
		end
	end

	# return
	return x_a, z_a

end

function consumer_colonisations(h, M_inc, x_p_old, z_p_old, i, j, neigh, x_a_old, z_a_old, alpha, ca, x_a, z_a)

	# get plant partners in current patch
	p_part_sp = get_plant_partners_in_current_patch(h, M_inc, x_p_old, i, j)

	# check neighbours
	# if at least one plant partner present in current patch
	if length(p_part_sp) > 0
		x_a, z_a = compute_consumer_colonisation_prob(z_p_old, i, j, p_part_sp, neigh, x_a_old, h, z_a_old, alpha, ca, x_a, z_a)
	end

	# return
	return x_a, z_a
end

function resource_extinctions_and_colonisations(n_p, x_p_old, i, j, M_inc, x_a_old, z_p_old, z_a_old, ep, x_p, z_p, neigh, cp, alpha)

	# resource extinctions and colonisations
    for h=1:n_p

    	# if species present -> extinctions
        if x_p_old[i,j,h] === 1
    		x_p, z_p = resource_extinctions(h, i, j, x_p, z_p, ep)
        continue
        end

        # if species absent -> colonisations
        if x_p_old[i,j,h] == 0
        	x_p, z_p = resource_colonisations(neigh, x_p_old, z_a_old, z_p_old, cp, i, j, h, x_p, z_p, M_inc, x_a_old, alpha)
        end

    end

	# return
	return x_p, z_p

end

function consumer_extinctions_and_colonisations(n_a, x_a_old, i, j, ea, x_a, z_a, M_inc, x_p_old, z_p_old, neigh, z_a_old, alpha, ca)

	# consumer extinctions and colonisations
	for h=1:n_a

    	# if species present -> extinctions
    	if x_a_old[i, j, h] === 1
    		x_a, z_a = consumer_extinctions(i, j, h, ea, x_a, z_a)
        	continue
    	end

    	# if species absent -> colonisations
    	if x_a_old[i,j,h] == 0
    	x_a, z_a = consumer_colonisations(h, M_inc, x_p_old, z_p_old, i, j, neigh, x_a_old, z_a_old, alpha, ca, x_a, z_a)
    	end
    end

	# return
	return x_a, z_a
end

function update_grid_after_colonisation_extinction(M_inc, p_current_sp, a_current_sp)

    # current incidence matrix
    current_inc = M_inc[p_current_sp, a_current_sp]

    # non-interacting species
    p_nonint_sp = Array{Int64}(undef, 0)
    a_nonint_sp = Array{Int64}(undef, 0)

    if sum(current_inc) === 0
    	p_nonint_sp = vcat(p_nonint_sp, p_current_sp)
        a_nonint_sp = vcat(a_nonint_sp, a_current_sp)
    end

	# return
	return current_inc, p_nonint_sp, a_nonint_sp

end
