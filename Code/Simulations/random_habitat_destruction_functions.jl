function habitat_destruction(n_patch, x_state, k, dn, n, x_p, x_a)

	# HABITAT DESTRUCTION
    # create matrix with patch numbers (col1) and states (col2)
    m_state = [1:n_patch  vec(x_state[:,:,k])]

    # matrix with pristine patches only
    m_pristine = m_state[m_state[:,2].!==0,:]

    # select dn patches at random and change their state to 0 (destroyed)
    patch_d = sample(m_pristine[:,1], dn, replace=false)
    m_state[in.(m_state[:,1], Ref(patch_d)),2] .= 0

    # new grid
    x_state[:,:,k+1] = reshape(m_state[:,2], (n,n))

    # species in destroyed patches become extinct
    x_p = x_p .* x_state[:,:,k+1]
    x_a = x_a .* x_state[:,:,k+1]

	# return
	return x_p, x_a, x_state


end
