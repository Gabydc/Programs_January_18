      function [fz]=koren_slope_limiter_fz(z,n,t,x,flux,dflux_dx)
%...
%...  The MatMol Group (2009)
%...
%...  function koren_slope_limiter returns the first derivative, fz, of a 
%...  flux function f(x) over the spatial domain z0 < z < zL from koren 
%...  slope limiter approximations.
%...
%...  The basic principle of slope limiters is to approximate the solution 
%...  x by combining a first-order, oscillation-free, monotone scheme with 
%...  a higher-order scheme in regions where the solution is sufficiently 
%...  smooth. In region of high-solution gradients, the first-order scheme 
%...  is prefered. 
%...
%...  The starting point of this procedure is the finite volume method in
%...  which the values x(i), defined on the local domain 
%...
%...    {z : z(i)< z <z(i+1)}, 
%...
%...  are the actual dependant variables ; x is supposed to be constant on 
%...  the local domain. To x(i) is associated the width of the local domain : 
%...
%...    dz(i) = (z(i+1) - z(i-1))/2.
%...
%...  Then, the PDE
%...
%...    x_t = - f_z(x,t)                                                                  (1)
%...
%...  is integrated by the method of lines following :
%...
%...    dx(i)     [f(x(i+1/2)) - f(x(i-1/2))]
%...    ----- = - ---------------------------
%...     dt                 dz(i)            
%...
%...  where
%...
%...    x(i+1/2) = x(z(i) + dz(i)/2)
%...
%...    x(i-1/2) = x(z(i) - dz(i)/2)
%...
%...  Because of his definition, x can take two values on  z(i)+dz(i)/2
%...  (equal to zR in the figure) and on z(i)-dz(i)/2 (equal to zL in 
%...  the figure) :
%...
%...           |           |                 |              |
%...           |           |                 |              |
%...           |           |                 |              |
%...           |           |                 |     x(i+1)   |
%...           |           |                 |--------------|
%...           |           |                 |              |
%...           |           |                 |              |
%...           |           |                 |              |
%...           |           |      x(i)       |              |
%...           |           |-----------------|              |
%...           |           |                 |              |
%...           |           |                 |              |
%...           |           |                 |              |
%...           |           |                 |              |
%...           |           |                 |              |
%...           |  x(i-1)   |                 |              |
%...           |-----------|                 |              |
%...        ___|___________|_________________|______________|_______
%...                      zL                zR
%...
%...  Moreover, the approximation of the solution to the piecewise constant
%...  values x(i) is only first order, and in order to construct a higher 
%...  (second) order approximation, the solution will be approximated by 
%...  piecewise linear functions. The slopes of those functions will be 
%...  limited by the slope limiter phi given in the koren case by
%...
%...    phi(r) = max(0, min(2*r, (1+2*r)/3, 2))
%...
%...  where r is the ratio of two consecutive solution derivatives : 
%...
%...  if the solution flows from left to right (f_x > 0), we have 
%...
%...             x(i+1)-x(i)
%...             -----------
%...             z(i+1)-z(i)
%...    r(i) = ---------------
%...             x(i)-x(i-1)
%...             -----------
%...             z(i)-z(i-1)
%...
%...  and, if the solution flows from right to left (f_x < 0),
%...
%...             x(i-1)-x(i)
%...             -----------
%...             z(i-1)-z(i)
%...    r(i) = ---------------
%...             x(i)-x(i+1)
%...             -----------
%...             z(i)-z(i+1)
%...
%...  It can be shown that x(i-1/2) and x(i+1/2) (called xL and xR in 
%...  the code) are given by
%...
%...  if f_x > 0 :
%...
%...                         phi(i-1)*dz(i-1)
%...    x(i-1/2) = x(i-1) +  ----------------*(x(i-1)-x(i-2))
%...                         dz(i-1)+dz(i-2)
%...
%...                       phi(i)*dz(i)
%...    x(i+1/2) = x(i) +  -------------*(x(i)-x(i-1))
%...                       dz(i)+dz(i-1)
%...
%...
%...  if f_x < 0 :
%...
%...                      phi(i)*dz(i)
%...    x(i-1/2) = x(i) - -------------*(x(i+1)-x(i))
%...                      dz(i+1)+dz(i)
%...
%...                        phi(i+1)*dz(i+1)
%...    x(i+1/2) = x(i+1) - ----------------*(x(i+2)-x(i+1))
%...                        dz(i+2)+dz(i+1)
%...
%...  argument list
%...
%...     z          independent variable (input)
%...
%...     n          number of grid points in the z domain including the
%...                boundary points (input)
%...
%...     t          time (input)
%...
%...     x          dependent variable (input)
%...
%...     flux       matlab function which computes f(x)
%...                call: f = fun(t,x) where where flux = 'fun'     (input)
%...
%...     dflux_dx   matlab function which computes the derivative of f with
%...                respect to x
%...                call: f_x = fun(t,x) where where dflux_dx = 'fun'h
%...                (in practice only the sign is used - there is no need 
%...                to compute an exact expression) (input)
%...
%...
      delta=1.0e-05;
      dz=zeros(n,1);
      dz(1)=(z(2)-z(1))/2;
      dz(2:n-1)=(z(3:n)-z(1:n-2))/2;
      dz(n)=(z(n)-z(n-1))/2;
      valdfdx=feval(dflux_dx,t,x);
%...
%...  computation of the fz derivative at the first left boundary point
      fz(1)=(feval(flux,t,x(2))-feval(flux,t,x(1)))/(z(2)-z(1));
%...
%...  computation of the fz derivative at the second left boundary point :
%...  fz(2) depends on the sign of valdfdx(2)
      if valdfdx(2) >= 0
          fz(2)=(feval(flux,t,x(2))-feval(flux,t,x(1)))/(z(2)-z(1));
      else
%...
%...      computation of r(2) and phi(2)
          if abs(x(2)-x(3)) < delta
              phi(2)=0;
          else
              r(2)=((x(1)-x(2))/(z(1)-z(2)))/((x(2)-x(3))/(z(2)-z(3)));
              phi(2)=max(0, min([2*r(2) (1+2*r(2))/3 2]));
          end
%...
%...      computation of r(3) and phi(3)
          if abs(x(3)-x(4)) < delta
              phi(3)=0;
          else
              r(3)=((x(2)-x(3))/(z(2)-z(3)))/((x(3)-x(4))/(z(3)-z(4)));
              phi(3)=max(0, min([2*r(3) (1+2*r(3))/3 2]));
          end
%...
%...      computation of xL and xR for the second left boundary point
          xL=x(2)-dz(2)*(x(3)-x(2))*phi(2)/(dz(2)+dz(3));
          xR=x(3)-dz(3)*(x(4)-x(3))*phi(3)/(dz(3)+dz(4));
          fz(2)=(feval(flux,t,xR)-feval(flux,t,xL))/dz(2);
      end
%...
%...  computation of the fz derivative at the interior points 3, ..., n-2 :
%...  fz(i) depends on the sign of valdfdx(i)
      for i=3:n-2
          if valdfdx(i) >= 0
%...
%...          computation of r(i) and phi(i)
              if abs(x(i)-x(i-1)) < delta
                  phi(i)=0;
              else
                  r(i)=((x(i+1)-x(i))/(z(i+1)-z(i)))/((x(i)-x(i-1))/(z(i)-z(i-1)));
                  phi(i)=max(0, min([2*r(i) (1+2*r(i))/3 2]));
              end
%...
%...          computation of r(i-1) and phi(i-1)
              if abs(x(i-1)-x(i-2)) < delta
                  phi(i-1)=0;
              else
                  r(i-1)=((x(i)-x(i-1))/(z(i)-z(i-1)))/((x(i-1)-x(i-2))/(z(i-1)-z(i-2)));
                  phi(i-1)=max(0, min([2*r(i-1) (1+2*r(i-1))/3 2]));
              end
%...
%...          computation of xL and xR for the interior points
              xL=x(i-1)+dz(i-1)*(x(i-1)-x(i-2))*phi(i-1)/(dz(i-1)+dz(i-2));
              xR=x(i)+dz(i)*(x(i)-x(i-1))*phi(i)/(dz(i)+dz(i-1));
              fz(i)=(feval(flux,t,xR)-feval(flux,t,xL))/dz(i);
          else
%...
%...          computation of r(i) and phi(i)
              if abs(x(i)-x(i+1)) < delta
                  phi(i)=0;
              else
                  r(i)=((x(i-1)-x(i))/(z(i-1)-z(i)))/((x(i)-x(i+1))/(z(i)-z(i+1)));
                  phi(i)=max(0, min([2*r(i) (1+2*r(i))/3 2]));
              end
%...
%...          computation of r(i+1) and phi(i+1)
              if abs(x(i+1)-x(i+2)) < delta
                  phi(i+1)=0;
              else
                  r(i+1)=((x(i)-x(i+1))/(z(i)-z(i+1)))/((x(i+1)-x(i+2))/(z(i+1)-z(i+2)));
                  phi(i+1)=max(0, min([2*r(i+1) (1+2*r(i+1))/3 2]));
              end
%...
%...          computation of xL and xR for the interior points
              xL=x(i)-dz(i)*(x(i+1)-x(i))*phi(i)/(dz(i)+dz(i+1));
              xR=x(i+1)-dz(i+1)*(x(i+2)-x(i+1))*phi(i+1)/(dz(i+1)+dz(i+2));
              fz(i)=(feval(flux,t,xR)-feval(flux,t,xL))/dz(i);
          end
      end
%...
%...  computation of the fz derivative at the penultimate point :
%...  fz(n-1) depends on the sign of valdfdx(n-1)
      if valdfdx(n-1) < 0
          fz(n-1)=(feval(flux,t,x(n))-feval(flux,t,x(n-1)))/(z(n)-z(n-1));
      else
%...
%...      computation of r(n-1) and phi(n-1)
          if abs(x(n-1)-x(n-2)) < delta
              phi(n-1)=0;
          else
              r(n-1)=((x(n)-x(n-1))/(z(n)-z(n-1)))/((x(n-1)-x(n-2))/(z(n-1)-z(n-2)));
              phi(n-1)=max(0, min([2*r(n-1) (1+2*r(n-1))/3 2]));
          end
%...
%...      computation of r(n-2) and phi(n-2)
          if abs(x(n-2)-x(n-3)) < delta
              phi(n-2)=0;
          else
              r(n-2)=((x(n-1)-x(n-2))/(z(n-1)-z(n-2)))/((x(n-2)-x(n-3))/(z(n-2)-z(n-3)));
              phi(n-2)=max(0, min([2*r(n-2) (1+2*r(n-2))/3 2]));
          end
%...
%...      computation of xL and xR for the penultimate point
          xL=x(n-2)+dz(n-2)*(x(n-2)-x(n-3))*phi(n-2)/(dz(n-2)+dz(n-3));
          xR=x(n-1)+dz(n-1)*(x(n-1)-x(n-2))*phi(n-1)/(dz(n-1)+dz(n-2));
          fz(n-1)=(feval(flux,t,xR)-feval(flux,t,xL))/dz(n-1);
      end
%...
%...  computation of the fz derivative at the right boundary point
      fz(n)=(feval(flux,t,x(n))-feval(flux,t,x(n-1)))/(z(n)-z(n-1));
      fz=fz';
