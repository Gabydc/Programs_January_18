%... The MatMol Group (2016)
     function Jac = jacobian_num2(t,x)
%...
%... numerical Jacobian matrix
     xt = bioreactor_pdes2(t,x);
     fac = [];
     thresh = 1e-6;
     threshv = thresh*ones(length(x),1);
     vectorized = 0;
     [Jac, fac] = numjac(@bioreactor_pdes2,t,x,xt,threshv,fac,vectorized);