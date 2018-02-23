%... The MatMol Group (2016)
     function Jac = jacobian_num(t,x)
%...
%... Set global variables
     global fac thresh vectorized
%...
%... numerical Jacobian matrix
     tstar = t;
     xstar = x;
     xtstar = fish_odes(tstar,xstar);
     threshv = [thresh;thresh;thresh];
     vectorized = 1;
     [Jac, fac] = numjac(@fish_odes,tstar,xstar,xtstar,threshv,fac,vectorized);
