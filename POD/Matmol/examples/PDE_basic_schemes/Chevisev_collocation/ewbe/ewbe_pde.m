%... The MatMol Group (2016)
     function xt = ewbe_pde(t,x,flag)

%... set global variables
     global a p mu delta c kappa phi0 z0;
     global zL zR z n D1;
%...
%... spatial derivatives
%...    
     xz = D1*x;
     xzz = D1*xz;
%...
%... temporal derivatives
%...
     xt = - a*(x.^p).*xz + delta*xzz;
     xt(1) = ewbe_exact(z(1),t) - x(1);
     xt(n) = ewbe_exact(z(n),t) - x(n);
