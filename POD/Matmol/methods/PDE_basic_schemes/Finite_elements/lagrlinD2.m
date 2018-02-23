%...  The Matmol group (2016)
    function out = lagrlinD2(h,n)
%...    
%... Computation of the second-order differentiation matrix of the 
%... finite element method with linear lagrangian elements

%... Main diagonal of the second-order differentiation matrix
    d0 = [0 -2*ones(1,n-2) 0];
%... Upper first diagonal of the second-order differentiation matrix
    dp1 = [0 ones(1,n-2)];
%... Lower first diagonal of the second-order differentiation matrix
    dm1 = [ones(1,n-2) 0];
%... Second order differentiation matrix
    out =sparse((1/h)*(diag(d0,0) + diag(dp1,1) + diag(dm1,-1)));
