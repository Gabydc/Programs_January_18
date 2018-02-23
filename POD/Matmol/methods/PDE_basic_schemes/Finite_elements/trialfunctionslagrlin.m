%...  The Matmol group (2016)
     function [phi_1 phi_2] = trialfunctionslagrlin(xi)
%...    
%... Computation of the lagrangian basis functions of the finite
%... element method
    phi_1 = (1-xi)/2;
    phi_2 = (1+xi)/2;
