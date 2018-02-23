%...  The Matmol group (2016)
     function out = hermitemass(h,n)
%...    
%... Computation of the mass matrix of the 
%... finite element method with quadratic hermitian elements

%... The number of variables is the double of discretization points
    N = 2*n;

%... Main diagonal of the mass matrix
    d0          = [78, 8, zeros(1,N-4), 78, 8];
    d0(3:2:N-3) = 156;
    d0(4:2:N-2) = 16;
    D0          = diag(d0,0);
%... Upper and lower first diagonals of the mass matrix
    d1          = [22, zeros(1,N-3), -22];
    d1(2:2:N-2) = 13;
    d1(3:2:N-3) = 0;
    Dp1          = diag(d1,1);
    Dm1          = diag(d1,-1);
%... Upper and lower second diagonals of the mass matrix
    d2          = zeros(1,N-2);
    d2(1:2:N-3) = 27;
    d2(2:2:N-2) = -6;
    Dp2          = diag(d2,2);
    Dm2          = diag(d2,-2);
%... Upper and lower third diagonals of the mass matrix
    d3          = zeros(1,N-3);
    d3(1:2:N-3) = -13;
    d3(2:2:N-4) = 0;
    Dp3          = diag(d3,3);
    Dm3          = diag(d3,-3);

%... Mass matrix
    out =sparse((h/210)*(D0+Dp1+Dm1+Dp2+Dm2+Dp3+Dm3));