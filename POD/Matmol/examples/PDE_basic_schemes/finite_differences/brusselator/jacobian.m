%... The MatMol Group (2016)
     function [Jac] = jacobian(t,x)
%...
%... global variables
     global A B alpha n dz
%...
%... Jacobian
     k(1:n-1) = 2*x(1:n-1).*x(1+n-1:n-1+n-1)-(B+1)-2*alpha;
     r(1:n-1) = B-2*x(1:n-1).*x(1+n-1:n-1+n-1);
     s(1:n-1) = x(1:n-1).^2;
     q(1:n-1) = -x(1:n-1).^2-2*alpha;
%...     
     d0 = [k q];
     dp1 = [alpha*ones(1,n-2) 0 alpha*ones(1,n-2)];
     dpnp2 = [s(1:n-1)];
     dm1 = dp1;
     dmnp2 = [r(1:n-1)];
%...     
     Jac = diag(d0)+diag(dp1,1)+diag(dm1,-1)+diag(dpnp2,n-1)+diag(dmnp2,-n+1);
