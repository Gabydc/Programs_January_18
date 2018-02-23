%... The MatMol Group (2016)

%... Korteweg-de Vries equation (solution with spectral methods)
%...
%... The Korteweg-de Vries equation is given by
%...
%...    u  = -a*u*u  - b*u                                  (1)
%...     t         z      zzz
%...
%... The propagation of a single soliton can be written as  
%...
%... u(z,t)= (3*c/a)*sech(0.5*sqrt(c/b)*(z-c*t))^2          (2)
%...
%... The initial conditions (ICs) for eq. (1) is derived
%... from (2)
%...
%... Dirichlet (BCs) for equations (1) are
%...
%...    u(0,t) = 0             u(L,t) = 0                   (3)
%...
%... The following code computes a solution to eqs. (1-3)
%...
     close all
     clear all
%...
%... Start a stopwatch timer
     tic
%...
%... Set global variables
     global a b c
     global z0 zL z n D1 D3
%...
%... Spatial grid
     z0 = -30.0;
     zL = 70.0;
     n = 256;               % for spectral method
%     n = 700;              % for FDs
%...
%... Initial, final times
     t0 = 0.0; tf = 25.0; 
%...  
%... Differentiation matrix
     [D1,z] = periodic_spectral_D1(z0,zL,n);
%... possibly compared with classical FDs
%     D1 = three_point_centered_D1(z);
%     D3 = seven_point_centered_D3(z);
%...
%... Initial conditions
%... u(z,t) = (3*c/a)*(sech(0.5*sqrt(c/b)*(z-c*t))^2
     a = 1;
     b = 1;
     c = 1;
     fz = 0.5*sqrt(c/b)*z;
     u = (3*c/a)*(sech(fz)).^2;
%...
%... Call to ODE solver
%...
%     method = 'leap_frog'
%...
%     method = 'RK4'
%...
%     method = 'ode45'
%...
%      method = 'ros23p'
%...
     method = 'ode15s'
%...
     switch method    
%...
%... leap_frog
      case('leap_frog')
      Dt = 0.001; Dtplot = 5.0;
      fzm1 = 0.5*sqrt(c/b)*(z+c*Dt);
      um1 = (3*c/a)*(sech(fzm1)).^2;
      [tout, uout] = leap_frog_solver(@kdv3_pde,t0,tf,u,um1,Dt,Dtplot);
%...
%... RK4
      case('RK4')
      Dt = 0.005; Dtplot  =5.0;
      [tout,uout] = rk4_solver(@kdv3_pde,t0,tf,u,Dt,Dtplot);
%...
%... ode45
      case('ode45')
      Dt = 5;
      options = odeset('RelTol',1e-6,'AbsTol',1e-6);
      [tout,uout] = ode45(@kdv3_pde,[t0:Dt:tf],u,options);
%...
%... ros23p
       case('ros23p')
       hmin = 0.0001;
       nstepsmax = 1e5;
       abstol = 1e-3;
       reltol = 1e-3;
       Dtplot = 1.0;
       [tout, uout, eout] = ros23p_solver(@kdv3_pde,@Jacobian,t0,tf,u,hmin,nstepsmax,abstol,reltol,Dtplot);
%...
%... ode15s
      case('ode15s')
      Dt = 1;
      options = odeset('RelTol',1e-3,'AbsTol',1e-8);
 %     options = odeset(options,'JPattern',jpattern(n));
 %     options = odeset(options,'Jacobian',@jacobian);
%...
 %     figure(1);
 %     ind = [1 round(n/4) round(n/2) round(3*n/4) n];
 %     options = odeset(options,'OutputFcn',@odeplot,'OutputSel',[ind]);
      [tout,uout] = ode15s(@kdv3_pde,[t0:Dt:tf],u,options);
%...
      end
%...
%...
     figure(2)
     plot(z,uout);
     hold on
     xlabel('z');
     ylabel('u(z,t)');
     title('Korteweg-de Vries equation')
%...
     figure(3)
     mesh(z,tout,uout);
	 shading interp;
     light;
     view(10,70);
     grid off;
     axis([z0 zL t0 tf min(min(uout)) max(max(uout))]);
     xlabel('z');
     ylabel('t');
     zlabel('u(t,z)');
%...
%... Read the stopwatch timer
     tcpu=toc;

  
