      function [fz]=kurg_centred_slope_limiter_fz(ne,n,z,t,x,flux,dflux_dx)
%...
%...  The MatMol Group (2009)
%...
%...  function kurg_centred_slope_limiter returns the first derivative, fz,
%...  of a vector flux function f(x) over the spatial domain z0 < z < zL
%...  from kurg_centred slope limiter approximations.
%...
%...  The basic principle of slope limiters is to approximate the solution 
%...  x by combining a first-order, oscillation-free, monotone scheme with 
%...  a higher-order scheme in regions where the solution is sufficiently 
%...  smooth. In region of high-solution gradients, the first-order scheme 
%...  is prefered. 
%...
%...  argument list
%...
%...     ne         number of scalar equations (input)
%...
%...     n          number of grid points in the z domain including the
%...                boundary points (input)
%...
%...     z          independent variable (input) : z(n)
%...
%...     t          time (input)
%...
%...     x          dependent variable (input) : x(n,ne)
%...
%...     flux       matlab function which computes f(x) : f(n,ne)
%...                call: f = fun(t,x) where flux = 'fun'     (input)
%...
%...     dflux_dx   matlab function which computes the derivative of f with
%...                respect to x : dfdx(n,ne,ne)
%...                call: f_x = fun(t,x) where dflux_dx = 'fun'
%...
%...
      for i=1:ne
          xtmp(1:n) = x(1:n,i);
          
          xtmpz(1) = (xtmp(2)-xtmp(1))/(z(2)-z(1));

          xtmpz1(1:n-2) = 2*(xtmp(3:n)-xtmp(2:n-1))./(z(3:n)-z(2:n-1))';
          xtmpz2(1:n-2) = (xtmp(3:n)-xtmp(1:n-2))./(z(3:n)-z(1:n-2))';          
          xtmpz3(1:n-2) = 2*(xtmp(2:n-1)-xtmp(1:n-2))./(z(2:n-1)-z(1:n-2))';
          testsgn       = sign(xtmpz1)+sign(xtmpz2)+sign(xtmpz3);
          for j=1:n-2
              if (testsgn(j) == 3) || (testsgn(j) == -3)
                  xtmpz(j+1)=sign(testsgn(j))*min([abs(xtmpz1(j)) abs(xtmpz2(j)) abs(xtmpz3(j))]);
              else
                  xtmpz(j+1)=0;
              end
          end

          xtmpz(n)    = (xtmp(n)-xtmp(n-1))/(z(n)-z(n-1));
          
          xtp(1:n-1)  = xtmp(2:n)-xtmpz(2:n).*(z(2:n)-z(1:n-1))'/2;
          xtm(1:n-1)  = xtmp(1:n-1)+xtmpz(1:n-1).*(z(2:n)-z(1:n-1))'/2;
          
          xp(1:n-1,i) = xtp(1:n-1);
          xm(1:n-1,i) = xtm(1:n-1);
      end
      
      dfludxp = feval(dflux_dx,n-1,ne,xp);
      dfludxm = feval(dflux_dx,n-1,ne,xm);

      
      for i=1:n-1
          size(dfludxp);
          ajac(1:ne,1:ne) = dfludxp(i,1:ne,1:ne);
          bajac           = balance(ajac);
          vp(1:ne)        = eig(bajac);
          ajac(1:ne,1:ne) = dfludxm(i,1:ne,1:ne);
          bajac           = balance(ajac);
          vp(ne+1:2*ne)   = eig(bajac);
          aspeed(i,1)     = max(abs(vp));
      end

      flzp = feval(flux,n-1,ne,t,xp);
      flzm = feval(flux,n-1,ne,t,xm);
      
      for j=1:ne
          fz(1,j)=0;
          fz(2:n-1,j)=(flzp(2:n-1,j)-flzp(1:n-2,j)+flzm(2:n-1,j)-flzm(1:n-2,j)...
                       -aspeed(2:n-1,1).*(xp(2:n-1,j)-xm(2:n-1,j))...
                       +aspeed(1:n-2,1).*(xp(1:n-2,j)-xm(1:n-2,j)))./(z(3:n)-z(1:n-2));
          fz(n,j)=0;
      end
          