%... The MatMol Group (2016)
     function S = jpattern(n)

%... sparsity pattern of the Jacobian matrix (FDs)
%      e = ones(n,1);
%      S = spdiags([e e e e e e e], [-3 -2 -1 0 1 2 3], n, n);
%...
%... compute the Jacobian sparsity pattern numerically at (t_start, u_star)
%...
      t_star = 0;
      u_star = ones(n,1);
      ut_star = kdv3_pde(t_star,u_star);
      tresh = 1e-14*ones(n,1);
      fac = [];
      vectorized = 1;
 %...     
      [dfdu,fac] = numjac(@kdv3_pde,t_star,u_star,ut_star,tresh,fac,vectorized);
 %...     
 %... replace nonzero sparse matrix elements with ones
      dfdus = sparse(dfdu);
      S = spones(dfdus);
      figure(4);
      spy(dfdus);
      xlabel('dependent variables');
      ylabel('equations');
      [njac,mjac] = size(dfdu);
      ntotjac = njac*mjac;
      non_zero = nnz(dfdus);
      non_zero_percent = non_zero/ntotjac*100;
      stat = sprintf('Jacobian sparisity pattern - nonzeros %d (%.3f%%)',non_zero,non_zero_percent);
      title(stat);
      pause
