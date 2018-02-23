%...  The MatMol Group (2016)
     function [f] = flux(n,ne,t,x)
%...
%... Function....... : Computes the flux for the Buckley-Leverett example
%... Input...........: n:  number of discretisation points
%...                   ne: number of partial differential equations
%...                   t:  current time
%...                   x:  current differential values

     f = (x.^2)./(x.^2+((1-x).^2));
