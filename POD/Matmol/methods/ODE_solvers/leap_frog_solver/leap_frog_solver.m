     function [tout, xout] = leap_frog_solver(odefunction, t0, tf, x0, xm1, Dt, Dtplot)
%... The Matmol group (2016)
%...
%... leap_frog_solver.m solves first-order differential equations using the leap frog method.
%... [tout, xout] = leap_frog_solver('f',t0,tf,x0,xm1,Dt,Dtplot) integrates the system of differential 
%... equations xt = f(t,x) from t0 to tf with initial conditions x0 and xm1.  
%... 'f' is a string containing the name of an ODE file. Function f(t,x) must return a column vector. 
%... Each row in solution array yout corresponds to a value returned in column vector t.
%...   
%... Argument list
%...
%... f - string containing the name of the user-supplied problem
%...     call: xt = problem_name(t,x) where f = 'problem_name'
%...     t - independent variable (scalar)
%...     x - solution vector
%...     xt - returned derivative vector; xt(i) = dx(i)/dt
%... t0 - initial value of t
%... tf - final value of t
%... x0 - initial value vector (at t = t0)
%... xm1 - initial value vector (at t = t0-Dt)
%... Dt - time step size
%... Dtplot - plot interval
%... tout - returned integration points (column-vector).
%... xout - returned solution, one solution row-vector per tout-value
%...
%... Initialization
%...
     plotgap = round(Dtplot/Dt);         % number of computation steps within a plot step
     Dt = Dtplot/plotgap;                
     nplots = round((tf - t0)/Dtplot);   % number of plots
     t = t0;                             % initialise t
     xold = xm1;                         % initialise x
     x = x0;                             % initialise x
     tout = t0;                          % initialise output value
     xout = x0';                         % initialise output value
%...    
%... implement Leap-frog method
%...
     for i = 1:nplots,
         for j = 1:plotgap,
%...        
%... use MATLAB's feval function to access the function file
%...
             xnew = xold + 2*Dt*feval(odefunction, t, x);
             t = t + Dt;
             xold = x;
             x = xnew;
         end
%...
%... add latest result to the output arrays
%...
         tout = [tout ; t];
         xout = [xout ; x'];
%...   
     end
