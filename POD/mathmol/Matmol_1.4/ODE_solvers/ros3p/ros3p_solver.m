       function [tout, xout] = ros3p_solver(odefunction,jacobian,t0,tf,x0,Dt,Dtplot)
%...
%...  The MatMol Group (2009)
%...
%...   this function solves first-order differential equations using a third-order fixed-step Rosenbrock method
%...   [tout, xout] = rk4_solver(@f,@J,t0,tf,x0,Dt,Dtplot) 
%...   integrates the system of differential equations xt=f(t,x) from t0 
%...   to tf with initial conditions x0.  f is a string containing the name 
%...   of an ODE file.  Function f(t,x) must return a column vector.  
%...   'f' is a string containing the name of an ODE file.  Function f(t,x) must return a column vector.  
%...   'J' is a string containing the name of the function, which evaluates the Jacobian J=[df/dx]
%...
%...   Each row in solution array yout corresponds to a value returned in column vector t.
%...
%...   ros3p_solver.m solves first-order differential equations using a fixed-step implicit Rosenbrock
%...   method for a series of points along the solution by repeatedly calling function ssros3p for a 
%...   single Rosenbrock step. 
%...
%... Argument list
%...
%... f - String containing name of user-supplied problem description
%...       Call: xt = problem_name(t,x) where f = 'problem_name'
%...             t - independent variable (scalar)
%...             x - solution vector
%...             xt - returned derivative vector; xt(i) = dx(i)/dt
%...
%... J - String containing name of user-supplied Jacobian
%...       Call: Jac = Jacobian(t,x) where J = 'Jacobian'
%...             t      - independent variable (scalar)
%...             x      - Solution vector.
%...             Jac    - Returned Jacobian matrix; Jac(i,j) = df(i)/dx(j)
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
%... method parameters
	 gamma = 7.886751345948129e-01;
     a21 = 1.267949192431123;
     a31 = 1.267949192431123;
     a32 = 0.0;
     c21 =-1.607695154586736;
     c31 = -3.464101615137755;
     c32 = -1.732050807568877;
     alpha(1) = 0;
     alpha(2) = 1.0;
     alpha(3) = 1.0;
     m(1) = 2.000000000000000e+00;
     m(2) = 5.773502691896258e-01;
     m(3) = 4.226497308103742e-01;
%...    
%... step through the solution
%...
     for i = 1:nplots,
         for j = 1:plotgap,
%...
%...         Single ros3p step
             [t,x] = ssros3p(odefunction,jacobian,t0,x0,Dt,gamma,a21,a31,a32,c21,c31,c32,alpha,m);
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
%... Next Rosenbrock step
     end
