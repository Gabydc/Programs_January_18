%...  The MatMol Group (2016)
     function xt = buckley_leverett_pde(t,x)
%...
%... set global variables
     global eps;
     global ne z0 zL z n method D1;
%...
%... boundary conditions at z = 0
     x(1) = 1;
%...
%... boundary conditions at z = zL
     x(n) = 0;
%...
%... spatial derivatives
     switch method;    
%...
%... van Leer slope limiter approximation
       case('Van_Leer')
           fz=vanleer_slope_limiter_fz(n,z,t,x,'flux','dflux_dx');
%... superbee slope limiter approximation
       case('superbee')
           fz=superbee_slope_limiter_fz(n,z,t,x,'flux','dflux_dx');
%... smart slope limiter approximation
       case('smart')
           fz=smart_slope_limiter_fz(n,z,t,x,'flux','dflux_dx');
%... monotonized center slope limiter approximation
       case('mc')
           fz=mc_slope_limiter_fz(n,z,t,x,'flux','dflux_dx');
%... minmod slope limiter approximation
       case('minmod')
           fz=minmod_slope_limiter_fz(n,z,t,x,'flux','dflux_dx');
%... Koren slope limiter approximation
       case('Koren')
           fz=koren_slope_limiter_fz(n,z,t,x,'flux','dflux_dx');
%... Kurganov-Tadmor centered approximation
       case('KT')
           fz=KT_centered_limiter_fz(ne,n,z,t,x,'fluxKT','dfluxKT_dx');
%...
     end
%...
%... second-order spatial derivative - direct differentiation
     nu = 4*x.*(1-x);
     diffusion = eps*D1*(nu.*(D1*x));
%...
%... temporal derivatives
%...
     xt = -fz + diffusion;
     xt(1) = 0;
     xt(n) = 0;
     t
