%... The MatMol Group (2016)
%...
%... Spring-Mass system
%...
%... consider a mass attached to a spring and assume that the motion takes
%... place on a table with no friction. In addition, there is no air resistance
%... and no dissipation of energy in the spring.
%...
%... the equation describing the motion of the mass is given by
%...
%... m*z   + k*z = 0
%...    tt  
%...
%... where
%...
%... z      position of the mass on the table
%... m      mass
%... k      spring constant
%...
%... this equation can also be rewritten as a system of two first order
%... ODEs by introducing v = z
%...                          t
%...
%... z  = v
%...  t
%... v  = -k/m * z
%...  t
%...
%... the following code computes a solution to this problem
%...
     close all
     clear all
%...
%... start a stopwatch timer
     tic
%...
%... set global variables
     global m k
%...
%... model parameters
     m = 1;
     k = 0.1;
%...
%... initial conditions
     t0 = 0;
     tf = 100;
     z = 1;
     v = 0;
     x = [z v]';
%...  
%... call to ODE solver (comment/decomment one of the methods to select a solver)
%...
      method = 'leap-frog'           
%...
%      method = 'celf'           
%...
%      method = 'ode45'     
%...
%      method = 'ode15s'     
%...
     switch method    
%...
%... leap-frog
       case('leap-frog')
           Dt=0.01;
           Dtplot=0.1;
           xm1=x;
           [tout, yout] = leap_frog_solver(@spring_mass_odes,t0,tf,x,xm1,Dt,Dtplot);
%... circularly exact leap-frog
       case('celf')
           Dt=0.01; 
           Dtplot=0.1;
           [tout, yout] = celf_solver(@spring_mass_odes,t0,tf,x,Dt,Dtplot);
%... ode45
       case('ode45')
           options = odeset('RelTol',1e-3,'AbsTol',1e-3);
           Dtplot=0.1;
           t=[t0:Dtplot:tf];
           [tout, yout] = ode45(@spring_mass_odes,t,x,options);
%... ode15s
       case('ode15s')
           options = odeset('RelTol',1e-3,'AbsTol',1e-3);
           Dtplot=0.1;
           t=[t0:Dtplot:tf];
           [tout, yout] = ode15s(@spring_mass_odes,t,x,options);
%...
     end
%...
%... plot results
     figure(1)
     subplot(2,1,1)
     plot(tout,yout(:,1),'k');
     xlabel('t');
     ylabel('z(t)');
     title('Position');
     axis([0 90 -1 1]);
     subplot(2,1,2)
     plot(tout,yout(:,2),'k');
     xlabel('t');
     ylabel('v(t)');
     title('Velocity')
     axis([0 90 -0.35 0.35]);
     figure(2)
     plot(yout(:,1),yout(:,2),'k');
     xlabel('z(t)');
     ylabel('v(t)');
     title('Phase portrait')
     axis([-1.5 1.5 -0.4 0.4]);
%...
%... read the stopwatch timer
     tcpu=toc;

  
