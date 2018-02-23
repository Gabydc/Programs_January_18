% The MatMol Group (2016)...
%... equal-width-Burgers equation
%...
%...
%... x  + a*(x^p)*x  - delta*x   - mu*x    = 0
%...  t            z          zz       zzt
%...
%...
%... where
%...
%...  x       dependent variable
%...
%...  t       time
%...
%...  z       space
%...
%...  mu      positive real constant
%...
%...  p       1 for the classical EWBE
%...
%... Reference
%...
%... S. Hamdi, W.H. Enright, W.E. Schiesser, J.J. Gottlieb
%... Exact solutions and invariants of motion for general types of
%... regularized long wave equations
%... Mathematics and Computers in Simulation 65 (2004), 535-545
%...
%...
     close all
     clear all
%...
%... start a stopwatch timer
     tic
%...
%... set global variables
     global a p mu delta c kappa phi0 z0;
     global zL zR z n D1;
%...
%... spatial grid
     zL = -50.0;
     zR = 50.0;
     n = 300;       % spectral method
%     n = 500;      % FDs
%...
%... spectral differentiation
     [D1,z] = chebyshev_spectral_D1(zL,zR,n);
%...
%... possibly compared to FDs
%     D1 = three_point_centered_uni_D1(zL,zR,n);
%...
%... parameters and initial conditions
     p = 1;
     a = 1.5;
     mu = 2;
     delta = 0.05;
     z0 = 0;
     x = ewbe_exact(z,0);    
%...  
%... call to ODE solver 
%...
%... ode15s
     t0 = 0;
     dt = 5;
     tf = 50;
     t = [t0:dt:tf];
     options = odeset('RelTol',1e-6,'AbsTol',1e-8);
     M = mass;
     options = odeset(options,'Mass',M,'MassSingular','yes');
%     options = odeset(options,'JPattern',jpattern(n));
%...
%     figure(1);
%     colordef black;
%     ind = [1 round(n/4) round(n/2) round(3*n/4) n];
%     options = odeset(options,'OutputFcn',@odeplot,'OutputSel',[ind]);
%...     
     [tout, xout] = ode15s(@ewbe_pde,t,x,options);
%...
%... plot results
     set(0,'defaultaxesfontsize',12,'defaultaxeslinewidth',.7,'defaultlinelinewidth',1.5,'defaultpatchlinewidth',.7);
%...
     figure(2);
     plot(z,xout,'k');
     xlabel('z');
     ylabel('x(t,z)');
     hold on
%...     
     for k=1:length(tout),
         for i=1:n,
             xexact(i) = ewbe_exact(z(i),tout(k));
         end;
         plot(z,xexact,'r')
     end;
%...
     figure(3)
%     mesh(z,tout,xout);
     surf(z,tout,xout);
     view(10,70);
     grid off;
     axis([zL zR t0 tf min(min(xout)) max(max(xout))]);
     xlabel('z');
     ylabel('t');
     zlabel('x(t,z)');
%...
%... read the stopwatch timer
     tcpu=toc;

  
