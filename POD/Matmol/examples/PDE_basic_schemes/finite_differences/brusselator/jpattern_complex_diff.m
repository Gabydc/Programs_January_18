%... The MatMol Group (2016)
     function S = jpattern_complex_diff
%...
%... Set global variables
     global n
%...
%... sparsity pattern of the Jacobian matrix based on a
%... numerical evaluation with complex differentiation
     tstar = 0;
     cstar = linspace(0.1,1,n);
     Tstar = linspace(300,700,n);
     xstar(1:2:2*n-1,1) = cstar';
     xstar(2:2:2*n,1) = Tstar';
     m = length(xstar);
     Jac = zeros(m,m);
     h = m*eps;
     for k = 1:m
         x1 = xstar;
         x1(k) = x1(k)+h*i;
         xt = brusselator_pdes(tstar,x1);
         Jac(:,k)=imag(xt)/h;
     end
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
