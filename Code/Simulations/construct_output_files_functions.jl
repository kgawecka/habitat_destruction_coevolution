function build_p_data(n, x, t_sim, n_patch)

	# create dataframe
	df_p = DataFrame(t = fill(t_sim, n_patch*n),
                            patch_no = repeat(1:n_patch, outer=n),
                            species = repeat(1:n, inner=n_patch),
                            p = reshape(x, (n_patch*n)))

	# return
	return df_p

end

function build_z_data(n, z, t_sim, n_patch)

	df_z = DataFrame(t = fill(t_sim, n_patch*n),
                            patch_no = repeat(1:n_patch, outer=n),
                            species = repeat(1:n, inner=n_patch),
                            z = reshape(z, (n_patch*n)))

	df_z = dropmissing(df_z)

	# return
	return df_z

end


function calculate_global_abundance_and_traits_and_store(x, z, n_patch, n, df_ab_traits, t_sim)

	# calculate global abundance
	p_new = sum(reshape(x, (n_patch, n)), dims=1) ./ n_patch

	# calculate global trait means and sd
	z_reshaped = reshape(z, (n_patch, n))
	z_mean = Vector{Union{Missing, Float32}}(missing, n)
	z_sd = Vector{Union{Missing, Float32}}(missing, n)
	for i=1:n
	    z_mean[i] = mean(skipmissing(z_reshaped[:,i]))
	    z_sd[i] = std(skipmissing(z_reshaped[:,i]))
	end

	# store results
	df_ab_traits[df_ab_traits.t.===t_sim,"ab"] .= p_new[1,:]
	df_ab_traits[df_ab_traits.t.===t_sim,"z_mean"] .= z_mean
	df_ab_traits[df_ab_traits.t.===t_sim,"z_sd"] .= z_sd

	# return
	return df_ab_traits


end
