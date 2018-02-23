%...  The Matmol group (2016)
     function [phi1, phi2, phi3, phi4] = trialfunctionsher2n(xi)
%...    
%... Computation of the hermitian basis functions of the finite
%... element method
    phi1 = ((xi.^2-3).*xi+2)/4;
    phi2 = (((xi-1).*xi-1).*xi+1)/4;
    phi3 = ((-xi.^2+3).*xi+2)/4;
    phi4 = (((xi+1).*xi-1).*xi-1)/4;
