      function [D]=five_point_biased_upwind_uni_D1(z0,zL,n,v)
%...
%...  The MatMol Group (2009)
%...
%...  function five_point_biased_upwind_uni_D1 returns the differentiation matrix 
%...  for computing the first derivative, xz, of a variable x over the spatial domain
%...  z0 < z < zL from biased-upwind five-point, fourth-order finite difference 
%...  approximations (this function replaces dss020)
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
%...  this function is an application of fourth-order directional
%...  differencing in the numerical method of lines.  it is intended
%...  specifically for the analysis of convective systems modelled by
%...  first-order hyperbolic partial differential equations as dis-
%...  cussed in subroutine bumd1d1.  the coefficients of the finite
%...  difference approximations used herein are taken from Bickley, W.
%...  G., Formulae for numerical differentiation, The Mathematical
%...  Gazette, pp. 19-27, 1941, n = 4, m = 1, p = 0, 1, 2, 3, 4.  the
%...  implementation is the **five-point biased upwind formula** of
%...  M. B. Carver and H. W. Hinds, The method of lines and the
%...  advection equation, Simulation, vol. 31, no. 2, pp. 59-69,
%...  august, 1978
%...
      dz=(zL-z0)/(n-1);
      r4fdz=1/(12.*dz);
%...
%...     (1)  finite difference approximation for positive v     
              if v > 0
%...
%...             sparse discretization matrix      
%...
%...             interior points      
                 e = ones(n,1);
                 D = spdiags([-1*e +6*e -18*e +10*e +3*e], [-3:1], n, n);
%...
%...             boundary points      
                 D([1 2 3],1:5) = [-25 +48 -36 +16 -3; -3 -10 +18 -6 +1; +1 -8 0 +8 -1];
                 D(n,n-4:n) = [+3 -16 +36 -48 +25];
              end;
%...
%...     (2)  finite difference approximation for negative v
              if v < 0
%...
%...             sparse discretization matrix      
%...
%...             interior points      
                 e = ones(n,1);
                 D = spdiags([-3*e -10*e +18*e -6*e +1*e], [-1:3], n, n);
%...
%...             boundary points      
                 D(1,1:5) = [-25 +48 -36 +16 -3];
                 D([n-2 n-1 n],n-4:n) = [+1 -8  0 +8 -1; -1 +6 -18 +10 +3; +3 -16 +36 -48 +25];
              end;                
%...      
      D=r4fdz*D;

