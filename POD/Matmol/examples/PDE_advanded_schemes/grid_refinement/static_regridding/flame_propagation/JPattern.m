%... The MatMol Group (2016)
     function S = jpattern(n)
%... sparsity pattern of the Jacobian matrix
     e = ones(2*n,1);
     S = spdiags([e e e e e e e e e e e], [-5 -4 -3 -2 -1 0 1 2 3 4 5], 2*n, 2*n);
