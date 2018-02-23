%...  The Matmol group (2016)
%...
%... Generalized equal-width wave equation (explicit solution with leap frog)
%...
%...
%... x  + a*(x^p)*x  - mu*x    = 0
%...  t            z       zzt
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
%...  p       1 for the classical EWWE
%...
%... Reference
%...
%... S. Hamdi, W.H. Enright, W.E. Schiesser, J.J. Gottlieb
%... Exact solutions and invariants of motion for general types of
%... regularized long wave equations
%... Mathematics and Computers in Simulation 65 (2004), 535-545
%...
     close all
     clear all
%...
%... start a stopwatch timer
     tic
%...
%... set global variables
     global a p mu c zi;
     global z0 zL z n D1;
     global L U
%...
%... spatial grid
     z0 = -15.0;
     zL = 70.0;
     n = 256;
%...
%... spectral differentiation
     [D1,z] = periodic_spectral_D1(z0,zL,n);
%...
%... initial conditions
     p = 1;
     a = 1;
     mu = 1;
     c = 0.3;
     zi = 0;
     x = (((p+1)*(p+2)*c*(sech(0.5*p*sqrt(1/mu)*z)).^2)/(2*a)).^(1/p);    
%...
%... Assemble mass matrix
     M = diag(ones(n,1),0)-mu*D1*D1;
     [L,U] = lu(M);
%...  
%... call to ODE solver 
%...
%... leap_frog
     t0 = 0;
     Dtplot = 5;
     tf = 120;
     Dt = 0.01;
     xm1 = ewwe_exact(z,-Dt);
     [tout, xout] = leap_frog_solver('ewwe_pde2',t0,tf,x,xm1,Dt,Dtplot);
%...
     set(0,'defaultaxesfontsize',12,'defaultaxeslinewidth',.7,'defaultlinelinewidth',.8,'defaultpatchlinewidth',.7);
%...
     figure(1)
     plot(z,xout,'k');
     xlabel('z');
     ylabel('x(t,z)');
     hold on
%...     
     for k=1:length(tout),
         xexact = ewwe_exact(z,tout(k));
         plot(z,xexact,':r')
     end;
%...
     figure(2)
     mesh(z,tout,xout);
     view(10,70);
     grid off;
     axis([z0 zL t0 tf min(min(xout)) max(max(xout))]);
     xlabel('z');
     ylabel('t');
     zlabel('x(t,z)');
%...
%... read the stopwatch timer
     tcpu=toc;

  
