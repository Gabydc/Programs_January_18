%... The MatMol Group (2016)
%...
%... Logistic equation
%...
%... N  = (a-b*N)*N
%...  t  
%...
%... where
%...
%...  a             growth rate without environmental influences
%...
%...  b             parameter representing the effect of increased population density
%...
%...  t             time
%...
%... the following code computes a solution to this problem
%...
     clear all
     close all
%...
%... set global variables
     global a b K N0
%...
%... model parameters
     a = 1;
     b = 0.5e-4;
     K = a/b;
%...
%... initial conditions
     t0 = 0;
     N0 = 1000;
     tf = 15;
     Dt = 0.1;
     Dtplot = 0.5;
     abstol = 1e-3;
     reltol = 1e-3;
     hmin = 1e-3;
     nstepsmax = 1000;
%...  
%... call to ODE solver
%...
     options = odeset('RelTol',1e-3,'AbsTol',1e-3);
     t = [t0:Dtplot:tf];
%     [tout, xout] = euler_solver(@logistic_ode,t0,tf,N0,Dt,Dtplot);
%     [tout, xout] = midpoint_solver(@logistic_ode,t0,tf,N0,Dt,Dtplot);
%     [tout, xout] = heun_solver(@logistic_ode,t0,tf,N0,Dt,Dtplot);
%     [tout, xout] = heun12_solver(@logistic_ode,t0,tf,N0,hmin,nstepsmax,abstol,reltol,Dtplot);
%     [tout, xout] = rk4_solver(@logistic_ode,t0,tf,N0,Dt,Dtplot);
%     [tout, xout] = rkf45_solver(@logistic_ode,t0,tf,N0,hmin,nstepsmax,abstol,reltol,Dtplot);
%     [tout, xout] = ros3p_solver(@logistic_ode,@jacobian_x,t0,tf,N0,Dt,Dtplot);
     [tout, xout] = ros23p_solver(@logistic_ode,@jacobian_x,@jacobian_t,t0,tf,N0,hmin,nstepsmax,abstol,reltol,Dtplot);
%...
%... plot results
     figure(1)
     plot(tout,xout,'k');
     hold on
     Nexact = logistic_exact(tout);
     plot(tout,Nexact,':r')
     xlabel('t');
     ylabel('N(t)');
     title('Logistic equation');