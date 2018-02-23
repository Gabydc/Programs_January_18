%... The MatMol Group (2016)
    function xt = flame_adaptive_pdes(t,x)     
%...
%... Set global variables
     global z0 zL nz D1;
%...
%... Transfer dependent variables
     rho = x(1:2:2*nz-1);
     T = x(2:2:2*nz);
%...
%... Boundary term
     if t < 2e-4,
        ignition = 0.2 + t/2e-4;
     else
        ignition = 1.2;
     end
%...
%... Source term
     NDA = 3.52e6*exp(-4./T);
%...
%... Spatial derivatives and boundary conditions
     rhoz = D1*rho;
     rhoz(1) = 0;
     rhoz(nz) = 0;
     rhozz = D1*rhoz;
     Tz = D1*T;
     Tz(1) = 0;
     Tzz = D1*Tz;
%...
%... Partial differential equations
%...
     rhot = rhozz - NDA.*rho;
     Tt = Tzz + NDA.*rho;
     Tt(nz) = 1e8*(ignition - T(nz));
%...
%... Transfer temporal derivatives
%...
     xt(1:2:2*nz-1) = rhot;
     xt(2:2:2*nz) = Tt;
     xt = xt';
