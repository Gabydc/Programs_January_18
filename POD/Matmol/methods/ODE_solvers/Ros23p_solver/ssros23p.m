     function [t,x,e] = ssros23p(odefunction,jacobian,time_derivative,t0,x0,h,gamma,a21,a31,a32,c21,c31,c32,d,alpha,m,mc)
%... The Matmol group (2016)
%...
%... Function ssros3p computes an ODE solution by an implicit third-order
%... Rosenbrock method for one step along the solution (by calls to 'odefunction' 
%... to define the ODE derivative vector, calls to 'jacobian' to define the Jacobian
%... and calls to time_derivative if the problem is non autonomous).  
%...
%... Argument list
%...
%... odefunction - string containing name of user-supplied problem
%... jacobian - string containing name of user-supplied Jacobian
%... time_derivative - sting containing name of user-supplied function time
%...                   derivative
%...
%... t0 - initial value of independent variable
%... x0 - initial condition vector
%... h - integration step
%... t - independent variable (scalar)
%... x - solution vector after one rkf45 step
%... e - estimate of truncation error of the solution vector
%...
%...   gamma,a21,a31,a32,c21,c31,c32,alpha,d,m,mc are the method parameters
%...
%... Jacobian matrix at initial (base) point
     [Jac] = feval(jacobian, t0, x0);
%...     
%... Time derivative at initial (base) point
     [Ft] = feval(time_derivative, t0, x0);
%...
%... Build coefficient matrix and perform L-U decomposition
	 CM = diag(1/(gamma*h)*ones(length(x0),1)) - Jac;
     [L,U] = lu(CM);
%...
%... stage 1
     xs = x0;
     [xt] = feval(odefunction, t0+alpha(1)*h, xs);
     rhs = xt + h*d(1)*Ft;
     xk1 = U\(L\rhs);
%...
%... stage 2
     xs = x0 + a21*xk1;
     [xt] = feval(odefunction, t0+alpha(2)*h, xs);
     rhs = xt + (c21/h)*xk1 + h*d(2)*Ft;
     xk2 = U\(L\rhs);
%...
%... stage 3
     xs = x0 + a31*xk1 + a32*xk2;
     [xt] = feval(odefunction, t0+alpha(3)*h, xs);
     rhs = xt + (c31/h)*xk1 + (c32/h)*xk2 + h*d(3)*Ft;
     xk3 = U\(L\rhs);
%...
%... second-order step
     x2 = x0 + mc(1)*xk1 + mc(2)*xk2 + mc(3)*xk3;
%...
%... third-order step
     x3  = x0 + m(1)*xk1 + m(2)*xk2 + m(3)*xk3;
%...
%... error evaluation
     t = t0 + h;
     e = x3 - x2;
     x = x3;
