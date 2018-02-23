%... The MatMol Group (2016)
%... Advection equation (solution with spectral methods)
%...
%...
%... c   = -v*c  + D*c   
%...  t        z      zz
%...
%...
%... where
%...
%...  c       concentration
%...
%...  t       time
%...
%...  z       space
%...
%...  v       fluid velocity
%...
%...  D       diffusion coefficient
%...
%...
%... the boundary conditions (BCs) are taken as
%...
%... c(z=z0,t) = c(t)
%...
%... c(z=zL,t) = 0               
%...  z
%...
%... where
%...
%...  [0,zL]    spatial domain
%...
%... the initial condition (IC) for Burgers' equations are taken as
%...
%... c(z,t=0) = 0
%...
%... the following code computes a solution to the advection equation
     clear all

%... start a stopwatch timer
     tic

%... set global variables
     global v D;
     global z0 zL n D1;

%... spatial grid
     z0 = 0.0;
     zL = 1.0;
     n  = 200;

%... spectral differentiation
     [D1,z] = chebyshev_spectral_D1(z0,zL,n);

%... model parameters
     v = 2;
     % v  = 2.*z;
     D = 0;

%... initial conditions
     c = exp(-100*(z-0.1))./(1 + exp(-100*(z-0.1)));

%... call to ODE solver
     %t            = 0:1:9;
     t            = 0:0.1:0.8;
     options      = odeset('RelTol',1e-6,'AbsTol',1e-6);
     [tout, yout] = ode15s(@advection_pde,t,c,options);

%... Figures
     plot(z,yout,'b');
     hold on
     xlabel('z');
     ylabel('c(t,z)');
     hold on

%... read the stopwatch timer
     tcpu=toc;


