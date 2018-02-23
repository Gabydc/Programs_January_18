%... The MatMol Group (2016)
     function xt = burgers_pde(t,x)
     t
%...
%... set global variables
     global mu;
     global z0 zL n dz D1 D2;
%...
%... boundary conditions at z = 0
     x(1) = burgers_exact(z0,t);
%...
%... boundary conditions at z = zL
     x(n) = burgers_exact(zL,t);
%...
%... second-order spatial derivative
%...    
     xzz = D2*x;
%...
%... first-order spatial derivative
     f = 0.5*x.^2;
     fminus = 0.5*(f(2:n-1)+f(1:n-2));
     fplus  = 0.5*(f(2:n-1)+f(3:n));
     fz = (fplus-fminus)/dz;
%...
%... temporal derivatives
%...
     xt(2:n-1,1) = -fz + mu*xzz(2:end-1);
     xt(1,1) = 0;
     xt(n,1) = 0;
