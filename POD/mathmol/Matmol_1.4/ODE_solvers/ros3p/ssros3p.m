    function [t,x] = ssros3p(odefunction,jacobian,t0,x0,Dt,gamma,a21,a31,a32,c21,c31,c32,alpha,m)
%...
%...  The MatMol Group (2009)
%...
%... Function ssros3p computes an ODE solution by an implicit third-order
%... Rosenbrock method for one step along the solution (by calls to 'odefunction' 
%... to define the ODE derivative vector and calls to 'jacobian' to define the Jacobian).  
%...
%... Argument list
%...
%... odefunction    - String containing name of user-supplied problem description
%...                        Call: xdot = odefunction(t,x)
%...                        t      - independent variable (scalar)
%...                        x      - Solution vector.
%...                        xdot   - Returned derivative vector; xdot(i) = dx(i)/dt
%...
%... jacobian      - String containing name of user-supplied Jacobian
%...                        Call: Jac = Jacobian(t,x)
%...                        t      - independent variable (scalar)
%...                        x      - Solution vector.
%...                        Jac    - Returned Jacobian matrix; Jac(i,j) = df(i)/dx(j)
%...
%...   t0          - initial value of independent variable
%...   x0          - initial condition vector
%...   Dt          - integration step
%...
%...   t           - independent variable
%...   x           - ODE solution vector after one rkf45 step
%...
%...   gamma,a21,a31,a32,c21,c31,c32,alpha,m are the method parameters
%...
%... Jacobian matrix at initial (base) point
     [Jac] = feval(jacobian, t0, x0);
%...
%... Build coefficient matrix and perform L-U decomposition
	 CM = diag(1/(gamma*Dt)*ones(length(x0),1)) - Jac;
     [L,U] = lu(CM);
%...
%... stage 1
     xs = x0;
     [xt] = feval(odefunction, t0+alpha(1)*Dt, xs);
     rhs = xt;
     xn1 = U\(L\rhs);
%...
%... stage 2
     xs = x0 + a21*xn1;
     [xt] = feval(odefunction, t0+alpha(2)*Dt, xs);
     rhs = xt + (c21/Dt)*xn1;
     xn2 = U\(L\rhs);
%...
%... stage 3
     xs = x0 + a31*xn1 + a32*xn2;
     [xt] = feval(odefunction, t0+alpha(3)*Dt, xs);
     rhs = xt + (c31/Dt)*xn1 + (c32/Dt)*xn2;
     xn3 = U\(L\rhs);
%...
%... third-order step
     x  = x0 + m(1)*xn1 + m(2)*xn2 + m(3)*xn3;
     t = t0 + Dt;
