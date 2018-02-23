%...  The Matmol group (2016)
     function M = mass(h,n,ne)
%...
%... Computation of the mass matrix of the finite element method

%... Main diagonal of the mass matrix
    d0 = repmat([2 4*ones(1,n-2) 2],1,ne);

%... First diagonal of the mass matrix
    d1 = [repmat([ones(1,n-1) 0],1,ne-1) ones(1,n-1)];

%... Mass matrix
    M = sparse((h/6)*(diag(d0,0) + diag(d1,1) + diag(d1,-1)));
