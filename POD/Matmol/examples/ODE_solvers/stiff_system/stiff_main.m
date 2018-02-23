%... The MatMol Group (2016)
%...
%... Introduction to stiff systems
%...
%... consider a system of two ODEs
%...
%... x(1)  = -a*x(1) + b*x(2)
%...     t  
%...
%... x(2)  = b*x(1) - a*x(2)
%...     t  
%...
%... the solution is given by
%...
%... x(1) = -exp(L1*t) + exp(L2*t)
%...
%... x(2) = exp(L1*t) + exp(L2*t)
%...
%... with L1 = -a - b et L2 = -a + b
%...
%... the initial conditions (ICs) are taken as
%...
%... x1(t=0) = 0    x2(t=0) = 2          
%...
%...
%... the following code computes a solution to this problem
%...
%... a = 50.5         L1 = -100        L1/L2=100
%... b = 49.5         L2 = -  1        (NONSTIFF)
%...
%... a =  500.5       L1 = -1000       L1/L2=1000
%... b =  499.5       L2 = -   1       (MODERATELY STIFF)
%...
%... a = 500000.5     L1 = -1000000    L1/L2=1000000
%... b = 499999.5     L2 = -      1    (VERY STIFF)
%...
     close all
     clear all
%...
%... start a stopwatch timer
     tic
%...
%... set global variables
     global a b
%...
%... model parameters
     cas = 1
%...
     if cas == 1
     a=1.5;
     b=0.5;
     
     elseif cas == 2
     a=50.5;
     b=49.5;
%...
     elseif cas == 3
     a=500.5;
     b=499.5;
%...
     elseif cas == 4
     a=500000.5;
     b=499999.5;
%...
     end
%...
%... initial conditions
     x=[0 2]';
%...  
%... call to ODE solver (comment/decomment one of the methods to select a solver)
%...
%      method = 'ode45'     
%...
      method = 'ode15s'     
%...
%      method = 'rkf45'     
%...
%      method = 'ros3p'     
%...
%      method = 'ros23p'     
%...
     switch method
%...         
%... ode45
       case('ode45')
           options = odeset('RelTol',1e-3,'AbsTol',1e-3);
           t=logspace(-5,1,100);
           [tout, yout] = ode45(@stiff_odes,t,x,options);
%... ode15s
       case('ode15s')
           options = odeset('RelTol',1e-3,'AbsTol',1e-3);
           t=logspace(-5,1,100);
           [tout, yout] = ode15s(@stiff_odes,t,x,options);
%... rkf45
       case('rkf45')
           t0=0;
           tf=10;
           hmin=0.0001;
           nstepsmax=1e5;
           abstol=1e-3;
           reltol=1e-3;
           Dtplot=0.01;
           [tout, yout, eout] = rkf45_solver(@stiff_odes,t0,tf,x,hmin,nstepsmax,abstol,reltol,Dtplot);
%... ros3p
       case('ros3p')
           t0=0;
           tf=10;
           Dt=0.001;
           Dtplot=0.01;
           [tout, yout] = ros3p_solver(@stiff_odes,@jacobian_x,t0,tf,x,Dt,Dtplot);
%... ros23p
       case('ros23p')
           t0=0;
           tf=10;
           hmin=0.0001;
           nstepsmax=1e5;
           abstol=1e-3;
           reltol=1e-3;
           Dtplot=0.01;
           [tout, yout, eout] = ros23p_solver(@stiff_odes,@jacobian_x,@jacobian_t,t0,tf,x,hmin,nstepsmax,abstol,reltol,Dtplot);
%...
     end
%...
%... plot results
     figure(1)
     plot(tout,yout(:,1),':r');
     xlabel('t');
     ylabel('x(t)');
     hold on
     plot(tout,yout(:,2),'--b');
%...
%... read the stopwatch timer
     tcpu=toc;

