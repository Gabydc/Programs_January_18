     function [tout, xout] = rk4_solver(odefunction,t0,tf,x0,Dt,Dtplot)
%... The Matmol group (2016)
%...
%... this function solves first-order differential equations using the
%... fourth-order fixed-step Runga-Kutta method
%... [tout, xout] = rk4_solver(@f,t0,tf,x0,Dt,Dtplot) 
%... integrates the system of differential equations xt=f(t,x) from t0 
%... to tf with initial conditions x0.  f is a string containing the name 
%... of an ODE file.  Function f(t,x) must return a column vector.  
%... Each row in solution array xout corresponds to a value returned 
%... in column vector t.
%...
%... rk4_solver.m solves first-order differential equations using 
%... the variable step RK Fehlberg 45 method for a series of points 
%... along the solution by repeatedly calling function ssrkf45 for a 
%... single RK Fehlberg 45 step.  The truncation error is estimated 
%... along the solution to adjust the integration step according to a 
%... specified error tolerance.
%...
%... Argument list
%...
%... f - String containing name of user-supplied problem description
%...       Call: xt = problem_name(t,x) where f = 'problem_name'
%...             t - independent variable (scalar)
%...             x - solution vector
%...             xt - returned derivative vector; xt(i) = dx(i)/dt
%...
%... t0 - initial value of t
%... tf - final value of t
%... x0 - initial value vector
%... Dt - integration interval
%... Dtplot - plot interval
%... tout - returned integration points (column-vector)
%... xout - returned solution, one solution row-vector per tout-value
%...
%... Initialization
%...
     plotgap = round(Dtplot/Dt);         % number of computation steps within a plot step
     Dt = Dtplot/plotgap;                
     nplots = round((tf - t0)/Dtplot);   % number of plots
     t = t0;                             % initialise t
     x = x0;                             % initialise x
     tout = t0;                          % initialise output value
     xout = x0';                         % initialise output value
%...    
%... step through the solution
%...
     for i = 1:nplots,
         for j = 1:plotgap,
%...
%...         RK 4 step
             [t,x] = ssrk4(odefunction,t0,x0,Dt);
%...
%...         Continue integration from new base point
             x0 = x;
             t0 = t;
%...
         end
%...
%... add latest result to the output arrays
%...
        tout = [tout ; t];
        xout = [xout ; x'];
%...
%... Next RK 4 step
     end
