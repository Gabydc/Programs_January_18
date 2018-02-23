%...  The Matmol group (2016)
    function [N1p N2p N3p N4p] = der1trialfunctionsher_1(x)
%...    
%... Derivative of the trial functions with hermite collocation    
    N1p = 3*(x^2-1)/4;
    N2p = ((3*x-2)*x-1)/4;
    N3p = -3*(x^2-1)/4;
    N4p = ((3*x+2)*x-1)/4;









