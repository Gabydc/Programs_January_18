     function [tout, yout] = celf_solver(odefunction, t0, tf, x0, Dt, Dtplot)
%... The Matmol group (2016)
%...
%... celf_solver solves first-order differential equations using the circularly exact leap frog (CELF)
%... method proposed by McLeod and Sanz-Serna (1982), Sanz-Serna and Manoranjan (1983).
%... [tout, yout] = celf_solver('f',t0,tf,x0,Dt,Dtplot) integrates the system of differential equations 
%... y'=f(t,x) from t0 to tf with initial conditions x0.  
%... 'f' is a string containing the name of an ODE file. Function f(t,x) must return a column vector. 
%... Each row in solution array yout corresponds to a value returned in column vector t.
%...
%... Argument list
%...
%... f        - String containing name of user-supplied problem description
%...                    Call: xdot = fun(t,x) where f = 'fun'
%...                    t      - independent variable (scalar)
%...                    x      - Solution vector.
%...                    xdot   - Returned derivative vector; xdot(i) = dx(i)/dt
%...
%... t0       - Initial value of t
%... tf       - Final value of t
%... x0       - Initial value vector (at t = t0)
%... Dt       - Initial step size
%... Dtplot   - Plot interval
%...
%... tout     - Returned integration points (column-vector).
%... yout     - Returned solution, one solution row-vector per tout-value.
%...
%...
%... Start integration
     t=t0;
     tini=t0;
     tout = t0;                                 %... initialize output value
     yout = x0';                                %... initialize output value
     nplots = round((tf - t0)/Dtplot);          %... number of plots
%...
     x = x0 + Dt*feval(odefunction, t, x0);     %... initialize with Euler's method
     xold = x0;
     t = t + Dt;
%...
%... add result to the output arrays
%...
     tout = [tout ; t];
     yout = [yout ; x'];
%...
     for i = 1:nplots,
       tplot = tini + i*Dtplot;
%...
%... While independent variable is less than the final
%... value, continue the integration
       while t <= tplot
%...        
%... use MATLAB's feval function to access the function file
%...
            f = feval(odefunction, t, x);
            h = abs((x - xold)'*f/(f'*f));
            xnew = xold + 2*h*f;
            t = t + h;
            xold = x;
            x = xnew;
%...           
       end
%...
%... add latest result to the output arrays and 
%... continue for  
        tout = [tout ; t];
        yout = [yout ; x'];
     end
