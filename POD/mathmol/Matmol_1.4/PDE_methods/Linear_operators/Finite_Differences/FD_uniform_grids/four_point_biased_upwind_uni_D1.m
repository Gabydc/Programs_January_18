      function [D]=four_point_biased_upwind_uni_D1(z0,zL,n,v)
%...
%...  The MatMol Group (2009)
%...
%...  function four_point_biased_upwind_uni_D1 returns the differentiation matrix 
%...  for computing the first derivative, xz, of a variable x over the spatial domain
%...  z0 < z < zL from biased-upwind four-point, third-order finite difference 
%...  approximations (this function replaces dss018)
%...
%...  argument list
%...
%...     z0      lower boundary value of z (input)
%...
%...     zL      upper boundary value of z (input)
%...
%...     n       number of grid points in the z domain including the
%...             boundary points (input)
%...
%...     v       fluid velocity (positive from left to right - only the sign is used) (input)
%...
%...  origin of the approximation
%...
%...  this function is an application of third-order directional
%...  differencing in the numerical method of lines.  It is intended
%...  specifically for the analysis of convective systems modelled by
%...  first-order hyperbolic partial differential equations as dis-
%...  cussed in function budm1d1.  the coefficients of the finite
%...  difference approximations used herein are taken from Bickley, W.
%...  G., Formulae for numerical differentiation, The Mathematical
%...  Gazette, pp. 19-27, 1941, n = 3, m = 1, p = 0, 1, 2, 3.  The
%...  implementation is the **four-point biased upwind formula** of
%...  M. B. carver and H. W. Hinds, The method of lines and the
%...  advection equation, Simulation, vol. 31, no. 2, pp. 59-69,
%...  August, 1978
%...
%...  compute the spatial increment
      dz=(zL-z0)/(n-1);
      r3fdz=1/(6*dz);
%...
%...     (1)  finite difference approximation for positive v     
              if v > 0
%...
%...             sparse discretization matrix      
%...
%...             interior points      
                 e = ones(n,1);
                 D = spdiags([+1*e -6*e +3*e +2*e], [-2:1], n, n);
%...
%...             boundary points      
                 D([1 2],1:4) = [-11 +18 -9 +2; -2 -3 +6 -1];
                 D(n,n-3:n) = [-2 +9 -18 +11];
              end;
%...
%...     (2)  finite difference approximation for negative v
              if v < 0
%...
%...             sparse discretization matrix      
%...
%...             interior points      
                 e = ones(n,1);
                 D = spdiags([-2*e -3*e +6*e -1*e], [-1:2], n, n);
%...
%...             boundary points      
                 D(1,1:4) = [-11 +18 -9 +2];
                 D([n-1 n],n-3:n) = [+1 -6 +3 +2; -2 +9 -18 +11];
              end;                
%...      
      D=r3fdz*D;
