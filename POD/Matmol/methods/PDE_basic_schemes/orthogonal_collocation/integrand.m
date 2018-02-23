%...  The Matmol group (2016)
    function out = integrand(x)
%...

global u1 n dz

%... Basis function
     [N1 N2 N3 N4]     = trialfunctionsher_1(x);
     [N1p N2p N3p N4p] = der1trialfunctionsher_1(x);

     out(1:n-1,1) = (2/dz)*(N1*u1(1:2:2*n-3) + N2*u1(2:2:2*n-2) + ...
                            N3*u1(3:2:2*n-1) + N4*u1(4:2:2*n)).*...
                            (N1p*u1(1:2:2*n-3) + N2p*u1(2:2:2*n-2)+...
                            N3p*u1(3:2:2*n-1) + N4p*u1(4:2:2*n));
