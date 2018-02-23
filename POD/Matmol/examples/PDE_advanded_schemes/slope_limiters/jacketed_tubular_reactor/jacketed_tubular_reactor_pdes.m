%... The MatMol Group (2016)
    function xt = jacketed_tubular_reactor_pdes(t,x)
%...
%... set global variables
    global Diff1 Diff2 v alpha delta gamma beta u uL uR
    global D1 D2 n dz nl z ne method
    t
%...
%... transfer variables
    x1(1:n,1) = x(1:n,1);
    x2(1:n,1) = x(n+1:2*n,1);
%...
%... impose boundary conditions
    x1(1,1) = x1(2,1)/(1+dz*v/Diff1);
    x1(n,1) = 2*x1(n-1,1)-x1(n-2,1);
    x2(1,1) = x2(2,1)/(1+dz*v/Diff2);
    x2(n,1) = 2*x2(n-1,1)-x2(n-2,1);
%...
%... spatial derivatives
     switch method;    
%...
%... van Leer slope limiter approximation
       case('Van_Leer')
           x1z=vanleer_slope_limiter_fz(n,z,t,x1,'flux','dflux_dx');
           x2z=vanleer_slope_limiter_fz(n,z,t,x2,'flux','dflux_dx');
%... superbee slope limiter approximation
       case('superbee')
           x1z=superbee_slope_limiter_fz(n,z,t,x1,'flux','dflux_dx');
           x2z=superbee_slope_limiter_fz(n,z,t,x2,'flux','dflux_dx');
%... smart slope limiter approximation
       case('smart')
           x1z=smart_slope_limiter_fz(n,z,t,x1,'flux','dflux_dx');
           x2z=smart_slope_limiter_fz(n,z,t,x2,'flux','dflux_dx');
%... monotonized center slope limiter approximation
       case('mc')
           x1z=mc_slope_limiter_fz(n,z,t,x1,'flux','dflux_dx');
           x2z=mc_slope_limiter_fz(n,z,t,x2,'flux','dflux_dx');
%... minmod slope limiter approximation
       case('minmod')
           x1z=minmod_slope_limiter_fz(n,z,t,x1,'flux','dflux_dx');
           x2z=minmod_slope_limiter_fz(n,z,t,x2,'flux','dflux_dx');
%... Koren slope limiter approximation
       case('Koren')
           x1z=koren_slope_limiter_fz(n,z,t,x1,'flux','dflux_dx');
           x2z=koren_slope_limiter_fz(n,z,t,x2,'flux','dflux_dx');
%... Kurganov-Tadmor centered approximation
       case('KT')
           xKT = [x1 x2];
           xz=KT_centered_limiter_fz(ne,n,z,t,xKT,@fluxKT,@dfluxKT_dx);
           x1z = xz(:,1);
           x2z = xz(:,2);
%...
     end
%...
%... second-order spatial derivative - direct differentiation
    x1zz = D2*x1;
    x2zz = D2*x2;
%...
%... jacket temperature
    u(1:nl,1) = uL;
    u(nl+1:n,1) = uR;
%...
%... temporal derivatives
%...
    x1t = Diff1*x1zz - v*x1z + alpha*delta*(1-x2).*exp(gamma*x1./(1+x1)) + beta*(u-x1);
    x2t = Diff2*x2zz - v*x2z + alpha*(1-x2).*exp(gamma*x1./(1+x1));
    x1t(1,1) = 0;
    x1t(n,1) = 0;
    x2t(1,1) = 0;
    x2t(n,1) = 0;
    xt = [x1t;x2t];