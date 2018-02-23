%... The MatMol Group (2016)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% implement the derivative df(x)/dx %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%...
     function [dfdx] = dflux_dx(n,ne,x)
%...
     global gam
     
     dfdx(1:n,1,1)=0;
     dfdx(1:n,1,2)=1;
     dfdx(1:n,1,3)=0;
     dfdx(1:n,2,1)=0.5*(gam-3)*((x(1:n,2)./x(1:n,1)).^2);
     dfdx(1:n,2,2)=-(gam-3)*x(1:n,2)./x(1:n,1);
     dfdx(1:n,2,3)=gam-1;
     dfdx(1:n,3,1)=-gam*(x(1:n,3).*x(1:n,2))./(x(1:n,1).^2)...
         +(gam-1)*((x(1:n,2)./x(1:n,1)).^3);
     dfdx(1:n,3,2)=gam*x(1:n,3)./x(1:n,1)-3*(gam-1)*((x(1:n,2)./x(1:n,1)).^2)/2;
     dfdx(1:n,3,3)=gam*x(1:n,2)./x(1:n,1);