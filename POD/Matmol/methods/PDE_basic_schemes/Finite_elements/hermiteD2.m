%...  The Matmol group (2016)
    function out = hermiteD2(h,n)
%...    
%... Computation of the second order differentiation matrix of the 
%... finite element method with quadratic hermitian elements

%... The number of variables is the double of discretization points
    N = 2*n;

%... Main diagonal of the second order differentiation matrix
    d0          = [-18, -8, zeros(1,N-4), -18, -8];
    d0(3:2:N-3) = -36;
    d0(4:2:N-2) = -16;
    D0          = diag(d0,0);
%... Upper first diagonals of the second order differentiation matrix
    d1          = [-33, zeros(1,N-3), 33];
    d1(2:2:N-1) = 3;
    d1(3:2:N-2) = 0;
    Dp1          = diag(d1,1);
%... Lower first diagonals of the second order differentiation matrix
    d1          = [-3, zeros(1,N-3), 3];
    d1(2:2:N-1) = 3;
    d1(3:2:N-2) = 0;
    Dm1         = diag(d1,-1);
%... Upper and lower second diagonals of the second order
%... differentiation matrix 
    d2          = zeros(1,N-2);
    d2(1:2:N-3) = 18;
    d2(2:2:N-2) = 2;
    Dp2         = diag(d2,2);
    Dm2         = diag(d2,-2);
%... Upper and lower third diagonals of the second order
%... differentiation matrix 
    d3          = zeros(1,N-3);
    d3(1:2:N-3) = -3;
    d3(2:2:N-4) = 0;
    Dp3          = diag(d3,3);
    Dm3          = diag(d3,-3);

%... Second order differentiation matrix
out =sparse((1/(15*h))*(D0+Dp1+Dm1+Dp2+Dm2+Dp3+Dm3));