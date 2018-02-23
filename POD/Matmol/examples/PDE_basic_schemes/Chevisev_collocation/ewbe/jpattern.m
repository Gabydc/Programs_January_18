%... The MatMol Group (2016)
     function S = jpattern(n)

%... compute the Jacobian sparsity pattern numerically at (t_start, x_star)
%...
     t_star = 0;
     x_star = ones(n,1);
     xt_star = ewbe_pde(t_star,x_star);
     tresh = 1e-14*ones(n,1);
     fac = [];
     vectorized = 1;
%...     
     [dfdx,fac] = numjac(@ewbe_pde,t_star,x_star,xt_star,tresh,fac,vectorized);
%...     
%... replace nonzero sparse matrix elements with ones
     dfdxs = sparse(dfdx);
     S = spones(dfdxs);
     figure(3);
     spy(dfdxs);
     xlabel('dependent variables');
     ylabel('equations');
     [njac,mjac] = size(dfdx);
     ntotjac = njac*mjac;
     non_zero = nnz(dfdxs);
     non_zero_percent = non_zero/ntotjac*100;
     stat = sprintf('Jacobian sparisity pattern - nonzeros %d (%.3f%%)',non_zero,non_zero_percent);
     title(stat);
     pause
