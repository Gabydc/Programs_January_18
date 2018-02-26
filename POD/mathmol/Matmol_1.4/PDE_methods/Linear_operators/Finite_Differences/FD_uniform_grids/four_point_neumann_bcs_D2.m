      function [xzz_z0,xzz_zL] = four_point_neumann_bcs_D2(z0,zL,x,xz_z0,xz_zL,n0,nL)
%...
%...  The MatMol Group (2009)
%...
%...  function four_point_neumann_bcs_D2 computes a four-point (2nd-order) approximation 
%...  ofthe second derivative at the boundaries, taking the value of the first 
%...  derivative at the boundary (Neumann boundary condition) into account.
%...
%...  argument list
%...
%...     x       one-dimensional array of the dependent variable
%...             to be differentiated (input)
%...
%...     xz_z0   value of the first derivative at the left boundary 
%...             z = z0 (input). (ignored if nl different from 1)
%...
%...     xz_zL   value of the first derivative at the right boundary 
%...             z = zL (input). (ignored if nu different from 1)
%...
%...     n0      integer flag for Neumann boundary condition at
%...             z = z0 (input).  (uxl=ux(1) is used if nl = 1)
%...
%...     nL      integer flag for Neumann of boundary condition at
%...             z = zL (input).  (uxu=ux(n) is used if nu = 1)
%...
%...     xzz_z0  value of the second derivative at the left boundary 
%...             z = z0 (output)
%...
%...     xzz_zL  value of the second derivative at the right boundary 
%...             z = zL (output)
%...
%...  origin of the approximations
%...
%...  the following derivation is for a set of second-order, four-point
%...  approximations for a second derivative that can be used at the
%...  boundaries of a spatial domain.  these approximations have the
%...  features
%...
%...     (1)  only interior and boundary points are used (i.e., no
%...          fictitious points are used)
%...
%...     (2)  the normal derivative at the boundary is included as part
%...          of the approximation for the second derivative
%...
%...     (3)  approximations for the boundary conditions are not used.
%...
%...  the derivation is by professor Gilbert A. Stengle, Department of
%...  Mathematics, Lehigh University, Bethlehem, PA 18015, and was done
%...  on december 7, 1985.
%...
%...  for an approximation at the left boundary, involving the points
%...  i, i+1, i+2 and i+3, consider the following taylor series expan-
%...  sions
%...
%...                  xz(i)( dz)   xzz(i)( dz)**2   xzzz(i)( dz)**3
%...  x(i+1) = x(i) + ---------- + -------------- + --------------- +...
%...                       1             2                 6
%...
%...
%...                  xz(i)(2dz)   xzz(i)(2dz)**2   xzzz(i)(2dz)**3
%...  x(i+2) = x(i) + ---------- + -------------- + --------------- +...
%...                       1             2                 6
%...
%...  if we now form the following linear combination, involving con-
%...  stants a, b, c and d to be determined, and use the preceding two
%...  taylor series,
%...
%...     a*x(i) + b*xz(i) + c*x(i+1) + d*x(i+2)
%...
%...  we have
%...
%...     a*x(i) + b*xz(i) + c*x(i+1) + d*x(i+2) =
%...
%...     (a + b + c + d)*x(i) +
%...
%...     (b + dz*c + 2*dz*d)*xz(i) +
%...
%...     (c*(dz**2)/2 + d*((2*dz)**2)/2)*xzz(i) +
%...
%...     (c*(dz**3)/6 + d*((2*dz)**3)/6)*xzzz(i) + o(dz**4)
%...
%...  the third derivative, xzzz(i), can be dropped by taking
%...
%...     c = -8*d
%...
%...  the second derivative, xzz(i), can be retained by taking
%...
%...     (dz**2)(c/2 + 2*d) = 1
%...
%...  which, when combined with the preceding result gives
%...
%...     d = -1/(2*(dz**2))
%...
%...     c = 4/(dz**2)
%...
%...  the first derivative, xz(i), can be dropped by taking
%...
%...     b + dz*c + 2*dz*d = 0
%...
%...  or
%...
%...     b = -dz*c - 2*dz*d = -4/dz - 2*dz*(-1/(2*(dz**2))) = -3/dz
%...
%...  finally, x(i), can be dropped by taking
%...
%...     a = - c - d = 8*d - d = -7*d = -7/(2*(dz**2))
%...
%...  if we now solve for the derivative of interest, xzz(i),
%...
%...     xzz(i) = -7/(2(dz**2))*x(i) - 3/dz*xz(i)
%...
%...              + 8/(dz**2)*x(i+1) - 1/(2*(dz**2))x(i+2) + o(dz**2)
%...
%...       = (1/(2*(dz**2)))*(-x(i+2) + 8*x(i+1) - 7*x(i) - 6*dz*xz(i))
%...
%...         + o(dz**2)
%...
%...  which is the four-point, second-order approximation for the second
%...  derivative, xzz(i), inclxding the first derivative, xz(i).
%...
%...  four checks of this approximation can easily be made for x(i) =
%...  1, x(i) = z, x(i) = z**2 and x(i) = z**3
%...
%...     xzz(i) = (1/(2*(dz**2)))*(-1 + 8*1 - 7*1 - 6*dz*0) = 0
%...
%...     xzz(i) = (1/(2*(dz**2)))*(-(z + 2*dz) + 8*(z + dz)
%...
%...              -7*z - 6*dz*1) = 0
%...
%...     xzz(i) = (1/(2*(dz**2)))*(-(z + 2*dz)**2 + 8*(z + dz)**2
%...
%...            - 7*(z**2) - 6*dz*(2*z))
%...
%...             = (-  z**2 -  4*z*dz - 4*dz**2
%...
%...               + 8*z**2 + 16*z*dz + 8*dz**2
%...
%...               - 7*z**2 - 12*z*dz)/(2*(dz**2)) = 2
%...
%...     xzz(i) = (1/(2*(dz**2)))*(-(z + 2*dz)**3 + 8*(z + dz)**3
%...
%...            - 7*(z**3) - 6*dz*(3*z**2))
%...
%...            = (1/(2*(dz**2)))*(- z**3 - 6*dz*z**2 - 12*z*dz**2
%...
%...            - 8*dz**3 + 8*z**3 + 24*dz*z**2 + 24*z*dz**2 + 8*dz**3
%...
%...            - 7*z**3 - 18*dz*z**2)
%...
%...            = (1/(2*(dz**2)))*(12*z*dz**2) = 6*z
%...
%...  the preceding approximation for xzz(i) can be applied at the
%...  left boundary value of z by taking i = 1.  an approximation at
%...  the right boundary is obtained by taking dz = -dz and reversing
%...  the subscripts in the preceding approximation, with i = n
%...
%...     xzz(i)
%...
%...       = (1/(2*(dz**2)))*(-x(i-2) + 8*x(i-1) - 7*x(i) + 6*dz*xz(i))
%...
%...         + o(dz**2)
%...
%...
%...  grid spacing
      n=length(x);
      dz=(z0-zL)/(n-1);
%...
%...  calculate xzz at the left boundary, including xz
      if n0==1
      xzz_z0 = (-7*x(1)+8*x(2)-1*x(3))/(2*dz^2)-6*xz_z0/(2*dz);
      end
%...
%...  calculate xzz at the right boundary, including xz
      if nL==1
      xzz_zL =(-7*x(n)+8*x(n-1)-1*x(n-2))/(2*dz^2)+6*xz_zL/(2*dz);
      end