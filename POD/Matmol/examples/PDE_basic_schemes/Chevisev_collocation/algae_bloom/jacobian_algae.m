%... The MatMol Group (2016)
    function Jac = jacobian_algae(t,x)

%... Set global variables
    global mu r K;
    global z0 zL z n D1;
    global tresh fac vectorized

%... analytical Jacobian matrix
%...	 Jac = mu*D1*D1 + diag(r*(1-x/K));
%...
%... numerical Jacobian matrix
    xt         = algae_bloom_pde(t,x);
    [Jac, fac] = numjac(@algae_bloom_pde,t,x,xt,tresh,fac,vectorized);