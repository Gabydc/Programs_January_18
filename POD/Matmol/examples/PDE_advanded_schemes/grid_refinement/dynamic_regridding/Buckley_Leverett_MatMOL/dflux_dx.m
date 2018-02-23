%...  The MatMol Group (2016)
     function [fx] = dflux_dx(n,ne,x)
%...
%... Function....... : Computes the flux gradient for the Buckley-Leverett
%...                   example
%... Input...........: n:  number of discretisation points
%...                   ne: number of partial differential equations
%...                   x:  current differential values

     fx = 2*x.*(1-x)./((x.^2+((1-x).^2)).^2);    
