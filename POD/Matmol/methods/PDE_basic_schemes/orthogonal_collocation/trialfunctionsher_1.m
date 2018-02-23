%...  The Matmol group (2016)
    function [N1 N2 N3 N4] = trialfunctionsher_1(x)
%...    
%... Trial functions with hermite collocation
    N1 = ((x^2-3)*x+2)/4;
    N2 = (((x-1)*x-1)*x+1)/4;
    N3 = ((-x^2+3)*x+2)/4;
    N4 = (((x+1)*x-1)*x-1)/4;

