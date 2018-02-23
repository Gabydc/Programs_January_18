%... The MatMol Group (2016)
%... The well known Burgers equation:
%...
%... x  = -x*x  + mu*x   
%...  t       z       zz
%...
%...
%... where
%...
%...  x      dependent variable
%...
%...  t       time
%...
%...  z       space
%...
%...  mu      diffusion coefficient
%...
%...
%... an exact solution is given by
%...
%... x(z,t) = burgers_exact(z,t)
%...
%... the boundary conditions (BCs) for Burgers' equation are taken as
%...
%... x(z=z0,t) = burgers_exact(z0,t)     (Dirichlet BCs)
%...
%... x(z=zL,t) = burgers_exact(zL,t)               
%...
%... where
%...
%...  zL      spatial domain
%...
%... the initial condition (IC) for Burgers' equations are taken as
%...
%... x(z,t=0) = burgers_exact(z,0)
%...
%... the following code computes a solution to Burgers' equation

     close all
     clear all

%... start a stopwatch timer
     tic

%... set global variables
     global mu;
     global z0 zL n D1 D2;

%... spatial grid
     z0 = 0.0;
     zL = 1.0;
     n  = 201;
     dz = (zL-z0)/(n-1);
     z  = [z0:dz:zL]';

%... model parameter
     % mu = 0.000001;
     mu = 0.001;

%... initial conditions
     for i=1:n,
         x(i) = burgers_exact_2(z(i),0);
     end;

%... select finite difference (FD) approximation of the spatial
%... derivative
     % method = 'centered_3';          % three point centered approximation
     % method = 'centered_5';          % five point centered approximation
     % method = 'upwind 2'           % two point upwind approximation
      method = 'upwind 3'           % two point upwind approximation
     % method = 'biased_upwind'    % five point, biased upwind approximation

     switch method
        %... five point centered approximation
         case('centered_5')
             D1 = five_point_centered_D1(z);
         %... three point centered approximation
         case('centered_3')
             D1=three_point_centered_D1(z);
         %... two point upwind approximation
         case('upwind 2')
             v = 1;
             D1=two_point_upwind_uni_D1(z0,zL,n,v);
         %... three point upwind approximation
         case('upwind 3')
             v = 1;
             D1=three_point_upwind_D1(z,v);
         %... five point, biased upwind approximation
         case('biased_upwind')
             v = 1;
             D1=five_point_biased_upwind_D1(z,v);
     end

     %... differentiation matrix (diffusive term)
     D2=five_point_centered_D2(z);

     %... call to ODE solver
     t = [0:0.1:1];
     options = odeset('RelTol',1e-6,'AbsTol',1e-6);
     [tout, yout] = ode45(@burgers_pde,t,x,options);

%... Plot the solution     
     plot(z,yout,'k');
     xlabel('z');
     ylabel('x(z,t)');
     axis([0,1,0,1.2]);

     title('Burgers equation')
     for k = 1:length(t),
         hold on
         for i = 1:n
             yexact(k,i) = burgers_exact_2(z(i),t(k));
         end
         plot(z,yexact,'r')
     end;

     figure
     hold on
     nc = 20;
     for k = 1:nc,
          zc(k) = z0+(k-1)*(zL-z0)/(nc-1);
          characteristics = (z-zc(k))/burgers_exact(zc(k),0);
          plot(z,characteristics,'r');
          xlabel('z')
          ylabel('t')
          axis([0 1 0 1])
     end;

%... read the stopwatch timer
     tcpu=toc;


