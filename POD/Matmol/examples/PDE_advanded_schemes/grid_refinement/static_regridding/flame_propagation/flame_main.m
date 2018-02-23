%... The MatMol Group (2016)
%...
%...  Flame propagation problem
%...
%... The classical Dwyer's and Sanders' flame propagation problem:
%...
%...    rho  = rho   - NDA*rho                                 (1)
%...       t      zz    
%...
%...    T  = T   + NDA*rho                                     (2)
%...     t    zz    
%...
%... where NDA = 3.52e6*exp(-4/T) and 
%...
%... The boundary conditions are given by
%...
%... rho (0,t) = 0                rho (1,t) = 0                (3)
%...    z                            z
%...
%... T (0,t) = 0                  T(1,t) = f(t)                (4)
%...  z
%...
%... where f(t) = 0.2 + t/2e-4  for t < 2e-4 and f(t) = 1.2 otherwise
%...
%... The initial conditions are given by
%...
%... rho(z,0) = 1                 T(z,0) = 0.2                  (5)
%...
%... The following code computes a solution to eqs. (1-5)
%...
     close all
     clear all
%...
%... Start a stopwatch timer
     tic
%...
%... Set global variables
     global z0 zL nz D1
     global nsteps maxsteps
%...
%... Spatial grid
     z0 = 0.0;
     zL = 1.0;
     nz = 401;
     dz = (zL-z0)/(nz-1);
     z = [z0:dz:zL]';
%...
%... Initial conditions
     rho = ones(nz,1);
     T = 0.2*ones(nz,1);
     x(1:2:2*nz-1) = rho;
     x(2:2:2*nz) = T;
     x = x';
%...
%... differentiation matrix
     D1 = three_point_centered_D1(z);
%...  
%... call to ODE solver 
%...     
     t0 = 0;
     tf = 0.006;
%...     
     Dt = tf/50;
     tspan = [t0:Dt:tf];
%...
%... do the integration 
        options = odeset('RelTol',1e-3,'AbsTol',1e-3);
        options = odeset(options,'JPattern',JPattern(nz));
        [t,y] = ode23s(@flame_adaptive_pdes,tspan,x,options);
%...
%... plot results
        figure(1)
        plot(z,y(:,1:2:2*nz-1),'b');
        xlabel('z')
        ylabel('\rho(z,t)')
        axis([0 1 0 1.2])
%...        
        figure(2)
        plot(z,y(:,2:2:2*nz),'r');
        xlabel('z')
        ylabel('T(z,t)')
        axis([0 1 0 1.5])
%...
%... read the stopwatch timer
%...
     tcpu=toc;
     
