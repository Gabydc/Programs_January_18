%...  The Matmol group (2016)
     function [phi1p, phi2p, phi3p, phi4p] = trialfunctionsher2np(x)
%...    
%... Computation of the derivative of the hermitian basis functions
%... of the finite element method
    phi1p = 3*(x^2-1)/4;
    phi2p = ((3*x-2)*x-1)/4;
    phi3p = -3*(x^2-1)/4;
    phi4p = ((3*x+2)*x-1)/4;
