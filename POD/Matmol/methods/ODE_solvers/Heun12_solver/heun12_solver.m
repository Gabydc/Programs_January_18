     function [tout, xout] = heun12_solver(odefunction, t0, tf, x0, hmin, nstepsmax, abstol, reltol, Dtplot)
%... The Matmol group (2016)
%...
%... this function solves first-order differential equations using Heun's method
%... [tout, xout] = heun12_solver(@f,t0,tf,x0,hmin,nstepsmax,abstol,reltol,Dtplot) 
%... integrates the system of differential equations xt=f(t,x) from t0 
%... to tf with initial conditions x0.  f is a string containing the name 
%... of an ODE file.  Function f(t,x) must return a column vector.  
%... Each row in solution array xout corresponds to a value returned 
%... in column vector t.
%...
%... heun12_solver.m solves first-order differential equations using 
%... a variable step method for a series of points along the solution 
%... by repeatedly calling function ssheun for a single integration step.
%... The truncation error is estimated along the solution to adjust the integration 
%... step according to a specified error tolerance.
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
%... hmin - minimum allowable time step
%... nstepsmax - maximum number of steps
%... abstol - absolute error tolerance
%... reltol - relative error tolerance
%... Dtplot - plot interval
%... tout - returned integration points (column-vector)
%... xout - returned solution, one solution row-vector per tout-value
%...
%... Start integration
     t = t0;
     tini = t0;
     tout = t0;                          %... initialize output value
     xout = x0';                         %... initialize output value
     eout = zeros(size(xout));           %... initialize output value
     nsteps = 0;                         %... initialize step counter
     nplots = round((tf - t0)/Dtplot);   %... number of plots
%...
%... Initial integration step
     h = 10*hmin;
%...
%... Step through nplots output points
     for i = 1:nplots
%...
%...   Final (output) value of the independent variable
       tplot = tini+i*Dtplot;
%...
%...   While independent variable is less than the final
%...   value, continue the integration
       while t <= tplot*0.9999
%...
%...     If the next step along the solution will go past the
%...     final value of the independent variable, set the
%...     step to the remaining distance to the final value
         if t+h > tplot, h = tplot-t; end
%...
%...     Single Heun step
         [t,x,e] = ssheun(odefunction,t0,x0,h);
%...
%...     Check if any of the ODEs have violated the error criteria
         if max( abs(e) > (abs(x)*reltol + abstol) )
%...
%...           Error violation, so integration is not complete.
%...           Reduce integration step because of error violation
%...           and repeat integration from the base point.
%...           Set logic variable for rejected integration step.
               h = h/2;
%...
%...     If the current step is less than the minimum allowable step,
%...     set the step to the minimum allowable value 
              if h < hmin,  h = hmin; end
%...          
%...     If there is no error violation, check if there is enough
%...     "error margin" to increase the integration step
         elseif max( abs(e) > (abs(x)*reltol + abstol)/4 )
%...
%...           The integration step cannot be increased, so leave 
%...           it unchanged and continue the integration from the 
%...           new base point.
               x0 = x; t0 = t;
%...
%...     There is no error violation and enough "security margin",
         else
%...
%...     double integration step and continue integration from new base point
              h = 2*h; x0 = x; t0 = t;
         end %if
%...
%... Continue while and check total number of integration steps taken
           nsteps=nsteps+1;
           if(nsteps > nstepsmax)  
              fprintf(' \n nstepsmax exceeded; integration terminated\n');
              break;
           end    
        end %while
%...
%... add latest result to the output arrays and continue for loop
        tout = [tout ; t];
        xout = [xout ; x'];
        eout = [eout ; e'];
     end % for
%...
%... End of Heun_solver
