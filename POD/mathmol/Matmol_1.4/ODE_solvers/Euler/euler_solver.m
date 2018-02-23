     function [tout,xout] = euler_solver(odefunction,t0,tf,x0,Dt,Dtplot)
%...
%...  The MatMol Group (2009)
%...
%... this function solves first-order differential equations using 
%... Euler's method.
%... [tout,xout] = euler_solver(@f,t0,tf,x0,Dt,Dtplot)
%... integrates the system of differential equations xt=f(t,x) from
%... t0 to tf with initial conditions x0.  f is a string containing
%... the name of an ODE file.  Function f(t,x) must return a column
%... vector.  Each row in solution array x corresponds to a value
%... returned in column vector t.
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
%... x0 - initial value vector
%... Dt - time step size
%... Dtplot - plot interval
%... tout - returned integration points (column-vector).
%... xout - returned solution, one solution row-vector per tout-value
%...
%... Initialization
     plotgap = round(Dtplot/Dt);      %... number of computation
                                      %... steps within a plot interval
     Dt = Dtplot/plotgap;
     nplots = round((tf-t0)/Dtplot);  %... number of plots
     t = t0;                          %... initialize t
     x = x0;                          %... initialize x
     tout = t0;                       %... initialize output value
     xout = x0';                      %... initialize output value
%...    
%... Implement Euler's method
     for i = 1:nplots
       for j = 1:plotgap
%...       
%...     Use MATLAB's feval function to access the function file,
%...     then take Euler step
         xnew = x + feval(odefunction,t,x)*Dt;
         t = t + Dt;
         x = xnew;
       end
%...
%...   Add latest result to the output arrays
       tout=[tout;t];
       xout=[xout;x'];
%...  
     end
