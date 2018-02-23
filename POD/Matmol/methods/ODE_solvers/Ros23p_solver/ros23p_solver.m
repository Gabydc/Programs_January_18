     function [tout, xout, eout] = ros23p_solver(odefunction,jacobian,time_derivative,t0,tf,x0,hmin,nstepsmax,abstol,reltol,Dtplot)
%... The Matmol group (2016)
%...
%... [tout, yout] = ros23p_solver('f','J','Ft',t0,tf,x0,hmin,nstepsmax,abstol,reltol,Dtplot) integrates  
%... a non-autonomous system of differential equations y'=f(t,x) from t0 to tf with initial conditions x0.  
%...
%... Each row in solution array xout corresponds to a value returned in column vector tout.
%... Each row in estimated error array eout corresponds to a value returned in column vector tout.
%...
%... ros23p_solver.m solves first-order differential equations using a variable-step implicit Rosenbrock
%... method for a series of points along the solution by repeatedly calling function ssros23p for a 
%... single Rosenbrock step. 
%...
%... The truncation error is estimated along the solution to adjust the integration step according to 
%... a specified error tolerance.
%...
%... Argument list
%...
%... f         - String containing name of user-supplied problem description
%...                    Call: xdot = fun(t,x) where f = 'fun'
%...                    t      - independent variable (scalar)
%...                    x      - Solution vector.
%...                    xdot   - Returned derivative vector; xdot(i) = dx(i)/dt
%...
%... J         - String containing name of user-supplied Jacobian
%...                    Call: Jac = fun(t,x) where J = 'fun'
%...                    t      - independent variable (scalar)
%...                    x      - Solution vector.
%...                    Jac    - Returned Jacobian matrix; Jac(i,j) = df(i)/dx(j)
%...
%... Ft        - String containing name of user-supplied function time
%....            derivative
%...                    Call: Ft = fun(t,x) where Ft = 'fun'
%...                    t      - independent variable (scalar)
%...                    x      - Solution vector.
%...                    Ft     - Returned time derivative vector; Ft(i) = df(i)/dt
%...
%... t0        - Initial value of t
%... tf        - Final value of t
%... x0        - Initial value vector
%... hmin      - minimum allowable time step
%... nstepsmax - maximum number of steps
%... abstol    - absolute error tolerance
%... reltol    - relative error tolerance
%... Dtplot    - Plot interval
%...
%... tout      - Returned integration points (column-vector).
%... xout      - Returned solution, one solution row-vector per tout-value.
%...
%... Initial integration step
     h = 10*hmin;
%...
%... method parameters
	 gamma = 0.5+sqrt(3)/6;
     a21 = 1.267949192431123;
     a31 = 1.267949192431123;
     a32 = 0.0;
     c21 =-1.607695154586736;
     c31 = -3.464101615137755;
     c32 = -1.732050807567788;
     d(1) = 7.886751345948129e-01;
     d(2) = -2.113248654051871e-01;
     d(3) = -1.077350269189626e+00;
     alpha(1) = 0;
     alpha(2) = 1.0;
     alpha(3) = 1.0;
     m(1) = 2.000000000000000e+00;
     m(2) = 5.773502691896258e-01;
     m(3) = 4.226497308103742e-01;
     mc(1) = 2.113248654051871e+00;
     mc(2) = 1.000000000000000e+00;
     mc(3) = 4.226497308103742e-01;
%... Start integration
     t = t0;
     tini = t0;
     tout = t0;                          %... initialize output value
     xout = x0';                         %... initialize output value
     eout = zeros(size(xout));           %... initialize output value
     nsteps = 0;                         %... initialize step counter
     nplots = round((tf-t0)/Dtplot);     %... number of outputs
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
%...   While independent variable is less than the final value,
%...   continue the integration
       while t <= tplot*0.9999
%...
%...     If the next step along the solution will go past the
%...     final value of the independent variable, set the
%...     step to the remaining distance to the final value
         if t+h > tplot, h = tplot-t; end
%...
%...     Single ros23p step
         [t,x,e] = ssros23p(odefunction,jacobian,time_derivative,t0,x0,h,gamma,a21,a31,a32,c21,c31,c32,d,alpha,m,mc);
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
         elseif max( abs(e) > (abs(x)*reltol + abstol)/8 )
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
%... End of ros23p_solver
