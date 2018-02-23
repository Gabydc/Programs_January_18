%...  The MatMol Group (2016)
     function [f] = flux(t,x)
%...
     f = (x.^2)./((x.^2)+(1-x).^2);