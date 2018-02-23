    function [t,x,e] = ssros23p(odefunction,jacobian,t0,x0,h,gamma,a21,a31,a32,c21,c31,c32,alpha,m,mc)
%...
%...  The MatMol Group (2009)
%...
%... Function ssros3p computes an ODE solution by an implicit third-order
%... Rosenbrock method for one step along the solution (by calls to 'odefunction' 
%... to define the ODE derivative vector and calls to 'jacobian' to define the Jacobian).  
%...
%... Argument list
%...
%... odefunction - string containing name of user-supplied problem
%... jacobian - string containing name of user-supplied Jacobian
%... t0 - initial value of independent variable
%... x0 - initial condition vector
%... h - integration step
%... t - independent variable (scalar)
%... x - solution vector after one rkf45 step
%... e - estimate of truncation error of the solution vector
%...
%...   gamma,a21,a31,a32,c21,c31,c32,alpha,m,mc are the method parameters
%...
%... Jacobian matrix at initial (base) point
     [Jac] = feval(jacobian, t0, x0);
%...
%... Build coefficient matrix and perform L-U decomposition
	 CM = diag(1/(gamma*h)*ones(length(x0),1)) - Jac;
     [L,U] = lu(CM);
%...
%... stage 1
     xs = x0;
     [xt] = feval(odefunction, t0+alpha(1)*h, xs);
     rhs = xt;
     xn1 = U\(L\rhs);
%...
%... stage 2
     xs = x0 + a21*xn1;
     [xt] = feval(odefunction, t0+alpha(2)*h, xs);
     rhs = xt + (c21/h)*xn1;
     xn2 = U\(L\rhs);
%...
%... stage 3
     xs = x0 + a31*xn1 + a32*xn2;
     [xt] = feval(odefunction, t0+alpha(3)*h, xs);
     rhs = xt + (c31/h)*xn1 + (c32/h)*xn2;
     xn3 = U\(L\rhs);
%...
%... second-order step
     x2 = x0 + mc(1)*xn1 + mc(2)*xn2 + mc(3)*xn3;
%...
%... third-order step
     x3  = x0 + m(1)*xn1 + m(2)*xn2 + m(3)*xn3;
%...
%... error evaluation
     t = t0 + h;
     e = x3 - x2;
     x = x3;
