%... The MatMol Group (2016)
     function Jac = jacobian_num(t,x)
%...
%... numerical Jacobian matrix
     tstar = t;
     xstar = x;
     xtstar = brusselator_pdes(tstar,xstar);
     fac = [];
     thresh = 1e-12;
     threshv = thresh*ones(length(x),1);
     vectorized = 0;
     [Jac, fac] = numjac(@brusselator_pdes,tstar,xstar,xtstar,threshv,fac,vectorized);