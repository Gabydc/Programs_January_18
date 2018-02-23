%... The MatMol Group (2016)
    function S = jpattern_algae(n)

%... sparsity pattern of the Jacobian matrix
    e = ones(n,1);
    S = spdiags([e e e e e], [-2 -1 0 1 2], n, n);
