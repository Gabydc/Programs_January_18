%...  The Matmol group (2016)
     function M = lagrlinmass(h,n,ne)
%...
%... Computation of the mass matrix of the finite element method 
%... with linear Lagrangian elements

%... Main diagonal of the mass matrix
    d0 = diag(repmat(([2 repmat(4,1,n-2) 2]),1,ne),0);

%... First upper diagonal of the mass matrix
    d1 = diag([repmat([ones(1,n-1) 0],1,ne-1) ones(1,n-1)],1);

%... First lower diagonal of the mass matrix
    dm1 = diag([repmat([ones(1,n-1) 0],1,ne-1) ones(1,n-1)],-1);

%... Mass matrix
    M = sparse((h/6)*(d0 + d1 + dm1));
