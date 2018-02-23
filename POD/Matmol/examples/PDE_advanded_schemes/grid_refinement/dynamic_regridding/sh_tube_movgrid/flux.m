%... The MatMol Group (2016)
     function [f] = flux(n,ne,t,x)
%...
     global gam
%...     
     f(:,1)=x(:,2);
     f(:,2)=(gam-1)*x(:,3)-0.5*(gam-3)*(x(:,2).^2)./x(:,1);
     f(:,3)=(gam*x(:,3)-0.5*(gam-1)*(x(:,2).^2)./x(:,1)).*(x(:,2)./x(:,1));