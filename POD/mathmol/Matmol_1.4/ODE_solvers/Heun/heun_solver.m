     function [tout, xout] = heun_solver(odefunction, t0, tf, x0, Dt, Dtplot)
%...
%...  The MatMol Group (2009)
%...
%... this function solves first-order differential equations using 
%... Heun's method
%... [tout, xout] = heun_solver(@f,t0, tf, x0, Dt, Dtplot) 
%... integrates the system of differential equations xt = f(t,x) from t0 
%... to tf with initial conditions x0.  f is a string containing the name 
%... of an ODE file.  Function f(t,x) must return a column vector.  
%... Each row in solution array xout corresponds to a value returned 
%... in column vector t.
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
%... implement Heun's method
%...
     for i = 1:nplots,
         for j = 1:plotgap,
%...        
%... use MATLAB's feval function to access the function file
%...
%...        Euler predictor
            f1 = feval(odefunction, t, x)*Dt;
            xEul = x + f1;
%...        Corrector            
            t = t + Dt;
            f2 = feval(odefunction, t, xEul)*Dt;
            x = x + (f1 + f2)/2;
         end
%...    
%... add latest result to the output arrays
%...
         tout = [tout ; t];
         xout = [xout ; x'];
%...    
     end
