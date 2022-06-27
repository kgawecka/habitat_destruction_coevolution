function simulate_coevolution(f, z, alpha, m, n_p_current, n_a_current, phi, z_p, z_a, theta, i, j, p_int_sp, a_int_sp)

	# number of species
	n_sp_current = n_p_current + n_a_current

	# trait change due to coevolution
	z_dif = transpose(f.*z) - f.*z
	q = f.*(exp.(-alpha.*(z_dif.^2)))
	q_n = q./sum(q,dims = 2)
	q_m = q_n.* m
	sel_dif = q_m .* z_dif
	r_mut = phi .* sum(sel_dif,dims=2)
	r_env = phi .* (1 .- m) .* (theta .- z)
	z_new = z .+ r_mut .+ r_env
	z_diff = abs.(z .- z_new)

    # update z_p and z_a
    z_p[i,j,p_int_sp] .= z_new[1:length(p_int_sp)]
    z_a[i,j,a_int_sp] .= z_new[(length(p_int_sp)+1):end]

    # return
	return z_p, z_a, q_m

end
