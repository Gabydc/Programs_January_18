%... The MatMol Group (2016)
     function S = jpattern(n)
%... sparsity pattern of the Jacobian matrix
     e = ones(n,1);
     S = spdiags([e e e e e e e], [-3 -2 -1 0 1 2 3], n, n);
