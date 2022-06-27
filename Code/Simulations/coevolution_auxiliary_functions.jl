function get_current_incidence_matrix(x_p, x_a, M_inc)

    # species present in current patch
    p_current_sp = findall(x_p[i, j, :] .=== 1)
    a_current_sp = findall(x_a[i, j, :] .=== 1)

    # check if at least one species present
    if length(p_current_sp) === 0 && length(a_current_sp) === 0
        #continue
    end

    # current incidence matrix
    current_inc = M_inc[p_current_sp, a_current_sp]

    # return
    return p_current_sp, a_current_sp, M_inc
end


function simulate_evolution(
    plant_parameters,
    animal_parameters,
    p_nonint_sp,
    a_nonint_sp,
    z_p,
    z_a,
    i,
    j,
)

    # coevolution parameters
    phi = vcat(
        plant_parameters[2, p_nonint_sp],
        animal_parameters[2, a_nonint_sp],
    )
    theta = vcat(
        plant_parameters[3, p_nonint_sp],
        animal_parameters[3, a_nonint_sp],
    )

    # trait values
    z = vcat(z_p[i, j, p_nonint_sp], z_a[i, j, a_nonint_sp])

    # trait change due to environmental selesciton only
    r_env = phi .* (theta .- z)
    z_new = z .+ r_env
    z_diff = abs.(z .- z_new)

    # update z_p and z_a
    z_p[i, j, p_nonint_sp] .= z_new[1:length(p_nonint_sp)]
    z_a[i, j, a_nonint_sp] .= z_new[(length(p_nonint_sp)+1):end]

    # return
    return z_p, z_a
end

function check_model_convergence(x_p_old, n_patch, n_p, x_p, x_a_old, n_a, x_a)

    p_p_old = sum(reshape(x_p_old, (n_patch, n_p)), dims = 1) ./ n_patch
    p_p_new = sum(reshape(x_p, (n_patch, n_p)), dims = 1) ./ n_patch
    p_a_old = sum(reshape(x_a_old, (n_patch, n_a)), dims = 1) ./ n_patch
    p_a_new = sum(reshape(x_a, (n_patch, n_a)), dims = 1) ./ n_patch
    dp = abs.(p_p_old - p_p_new) ./ p_p_old
    da = abs.(p_a_old - p_a_new) ./ p_a_old
    dp[isnan.(dp)] .= 0
    da[isnan.(da)] .= 0

    # return
    return dp, da

end
