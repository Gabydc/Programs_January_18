%... The MatMol Group (2016)
     function xt = brusselator_pdes2(t,x)
%...
%... global variables
     global A B alpha n dz ul ur vl vr
%...
%... separate u and v from x
     u(2:n) = x(1:n-1);
     v(2:n) = x(n:2*n-2);
     u(1) = ul;
     u(n+1) = ur;
     v(1) = vl;
     v(n+1) = vr;
%...
%..  PDEs
     ut(1:n-1) = A + (u(2:n).^2).*v(2:n) - (B+1)*u(2:n) + alpha*(u(3:n+1)-2*u(2:n)+u(1:n-1))/(dz^2);
     vt(1:n-1) = B*u(2:n)-(u(2:n).^2).*v(2:n) + alpha*(v(3:n+1)-2*v(2:n)+v(1:n-1))/(dz^2);
%...
%... transfer temporal derivatives
     xt = [ut  vt]';
     