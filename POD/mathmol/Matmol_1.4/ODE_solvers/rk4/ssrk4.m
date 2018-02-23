    function [t,x] = ssrk4(odefunction,t0,x0,Dt)
%...
%...  The MatMol Group (2009)
%...
%... Function ssrk4 computes an ODE solution by the classical fourth
%... order RK method for one step along the solution (by calls to 'odefunction' 
%... to define the ODE derivative vector).  
%...
%... Argument list
%...
%... odefunction - string containing name of user-supplied problem
%... t0 - initial value of independent variable
%... x0 - initial condition vector
%... h - integration step
%... t - independent variable (scalar)
%... x - solution vector after one rkf45 step
%...
%... Derivative vector at initial (base) point
     [xt0] = feval(odefunction, t0, x0);
%...
%... k1, advance of dependent variable vector and
%... independent variable for calculation of k2
     k1 = Dt*xt0;
     x = x0 + 0.5*k1;
     t = t0 + 0.5*Dt;
%...
%... Derivative vector at new x, t
     [xt] = feval(odefunction, t, x);
%...
%... k2, advance of dependent variable vector and
%... independent variable for calculation of k3
     k2 = Dt*xt;
     x = x0 + 0.5*k2;
     t = t0 + 0.5*Dt;
%...
%... Derivative vector at new x, t
     [xt] = feval(odefunction, t, x);
%...
%... k3, advance of dependent variable vector and
%... independent variable for calculation of k4
     k3 = Dt*xt;
     x = x0 + k3;
     t = t0 + Dt;
%...
%... Derivative vector at new x, t
     [xt] = feval(odefunction, t, x);
%...
%... k4
     k4 = Dt*xt;
%...
%... Fourth order step
     x = x0 + (1.0/6.0)*(k1 + 2.0*k2 + 2.0*k3 + k4);
     t = t0 + Dt;
