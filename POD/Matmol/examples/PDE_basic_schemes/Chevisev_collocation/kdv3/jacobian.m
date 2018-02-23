%... The MatMol Group (2016)
     function Jac = jacobian(t,u)

%... numerical Jacobian matrix    
     tstar = t;
     ustar = u;
     utstar = kdv3_pde(tstar,ustar);
     fac = [];
     thresh = 1e-12;
     threshv = thresh*ones(length(u),1);
     vectorized = 1;
     [Jac, fac] = numjac(@kdv3_pde,tstar,ustar,utstar,threshv,fac,vectorized);