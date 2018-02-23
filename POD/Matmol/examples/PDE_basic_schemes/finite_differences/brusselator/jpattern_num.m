%... The MatMol Group (2016)
     function S = jpattern_num
%...
%... Set global variables
     global n
%...
%... sparsity pattern of the Jacobian matrix based on a
%... numerical evaluation
     tstar = 0;
     cstar = linspace(0.1,1,n);
     Tstar = linspace(300,700,n);
     xstar(1:2:2*n-1,1) = cstar';
     xstar(2:2:2*n,1) = Tstar';
     xtstar = brusselator_pdes(tstar,xstar);
     fac = [];
     thresh = 1e-16;
     vectorized = 0;
     [Jac, fac] = numjac(@brusselator_pdes,tstar,xstar,xtstar,thresh,fac,vectorized);
%...
%... replace nonzero elements by "1" 
%... (so as to create a "0-1" map of the Jacobian matrix) 
     S = spones(sparse(Jac));
%...
%... plot the map
     figure
     spy(S);
     xlabel('dependent variables');
     ylabel('semi-discrete equations');
%...
%... compute the percentage of non-zero elements
     [njac,mjac] = size(S);
     ntotjac = njac*mjac;
     non_zero = nnz(S);
     non_zero_percent = non_zero/ntotjac*100;
     stat = sprintf('Jacobian sparsity pattern - nonzeros %d (%.3f%%)',non_zero,non_zero_percent);
     title(stat);
