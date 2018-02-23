%... The MatMol Group (2016)
     function xt = brusselator_pdes(t,x)
%...
%... global variables
     global A B alpha n dz ul ur vl vr
%...
%... separate u and v from x
     u = x(1:n);
     v = x(n+1:2*n);
%...
%... boundary conditions
     u(1) = ul; ut(1) = 0;
     u(n) = ur; ut(n) = 0;
     v(1) = vl; vt(1) = 0;
     v(n) = vr; vt(n) = 0;
%...
%..  PDEs
     ut(2:n-1) = A + (u(2:n-1).^2).*v(2:n-1) - (B+1)*u(2:n-1) + alpha*(u(3:n)-2*u(2:n-1)+u(1:n-2))/(dz^2);
     vt(2:n-1) = B*u(2:n-1)-(u(2:n-1).^2).*v(2:n-1) + alpha*(v(3:n)-2*v(2:n-1)+v(1:n-2))/(dz^2);
%...
%... transfer temporal derivatives
     xt = [ut  vt]';
     