     function [t,x,e] = ssheun(odefunction,t0,x0,h)
%...
%...  The MatMol Group (2009)
%...
%... Function ssheun computes an ODE solution by the classical Heun method
%... for one step along the solution (by calls to 'odefunction' to 
%... define the ODE derivative vector).  It also estimates the truncation 
%... error of the solution, and applies this estimate as a correction to the 
%... solution vector.
%...   
%... Argument list
%...
%... odefunction - string containing name of user-supplied problem
%... t0 - initial value of independent variable
%... x0 - initial condition vector
%... h - integration step
%... t - independent variable (scalar)
%... x - solution vector after one integration step
%... e - estimate of truncation error of the solution vector
%...
%... Euler predictor     
     f1 = feval(odefunction, t0, x0)*h;
     xEul = x0 + f1;
%...
%....Corrector
     t = t0 + h;
     f2 = feval(odefunction, t, xEul)*h;
%...
%... x = x0 + (f1+f2)/2
%... or alternatively    
     e = (f2 - f1)/2;
     x = xEul + e;
