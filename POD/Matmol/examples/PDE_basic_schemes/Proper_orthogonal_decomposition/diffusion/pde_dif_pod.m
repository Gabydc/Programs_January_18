%... The MatMol Group (2016)    
    function m_t = pde_dif_pod(t, m, A, kappa)

%... PDE
    m_t = kappa*A*m;