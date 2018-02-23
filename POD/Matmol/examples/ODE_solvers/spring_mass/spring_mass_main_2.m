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
%...  
%... call to ODE solver
%...
%... leap-frog (direct solution of second order equation)
           Dt=0.05;
           Dtplot=0.1;
           zm1=z;
           [tout, yout] = leap_frog_2_solver(@spring_mass_ode_2, t0, tf, z, zm1, Dt, Dtplot);
%...
%... plot results
     figure(1)
     plot(tout,yout(:,1));
     hold on
     xlabel('t');
     ylabel('z(t)');
     title('Position');
%...
%... read the stopwatch timer
     tcpu=toc;

  
