%...  The MatMol Group (2016)
     function [fx] = dfluxKT_dx(ne,t,x)
%...    
    fx = 2*x.*(1-x)./((x.^2)+(1-x).^2).^2;
