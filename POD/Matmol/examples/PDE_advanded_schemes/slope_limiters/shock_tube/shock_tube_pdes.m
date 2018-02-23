%... The MatMol Group (2016)
    function xt = shock_tube_pdes(t,x)
    global n z D2 eps npdes
    t
%...
%... separate dependent variables
%...
     u = zeros(n,npdes);
     u(1:n,1) = x(1:n,1);
     u(1:n,2) = x(n+1:2*n,1);
     u(1:n,3) = x(2*n+1:3*n,1);
%...
%... B.Cs.
%...
     u(1,1) = u(2,1);
     u(n,1) = u(n-1,1);
     u(1,2) = 0;
     u(n,2) = 0;
     u(1,3) = u(2,3);
     u(n,3) = u(n-1,3);
%...
%... spatial derivatives
%...
%... convection
%...
     fz = KT_centered_limiter_fz(npdes,n,z,t,u,@flux,@dflux_dx);
%...
%... diffusion
%...
     uzz = D2*u;
%...
%... temporal derivatives
%...
     xt(1:n,1) = - fz(:,1) + eps*uzz(:,1);
     xt(n+1:2*n,1) = - fz(:,2) + eps*uzz(:,2);
     xt(2*n+1:3*n,1) = - fz(:,3) + eps*uzz(:,3);
     xt(1,1) = 0;
     xt(n,1) = 0;
     xt(n+1,1) = 0;
     xt(2*n,1) = 0;
     xt(2*n+1,1) = 0;
     xt(3*n,1) = 0;
