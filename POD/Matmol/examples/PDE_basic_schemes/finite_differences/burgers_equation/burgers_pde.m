%... The MatMol Group (2016)
     function xt = burgers_pde(t,x)
%...
%... set global variables
     global mu;
     global z0 zL n D1 D2;
%...
%... boundary conditions at z = 0
     x(1) = burgers_exact(z0,t);
%...
%... boundary conditions at z = zL
     x(n) = burgers_exact(zL,t);
%...
%... second-order spatial derivative   
     xzz = D2*x;
%...
%... first-order spatial derivative
     xz = D1*x;
%...
%... temporal derivatives
%...
     xt = -x.*xz + mu*xzz;
     xt(1) = 0;
     xt(n) = 0;
