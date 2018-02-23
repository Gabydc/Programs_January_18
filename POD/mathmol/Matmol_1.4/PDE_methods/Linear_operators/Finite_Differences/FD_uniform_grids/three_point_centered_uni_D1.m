      function [D] = three_point_centered_uni_D1(z0,zL,n)
%...
%...  The MatMol Group (2009)
%...
%...  function three_point_centered_uni_D1 returns the differentiation matrix 
%...  for computing the first derivative, xz, of a variable x over the spatial domain
%...  z0 < z < zL from classical three-point, second-order finite difference 
%...  approximations (this function replaces dss002)
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
%...  origin of the approximation
%...
%...                                       2
%...  x1  = (1/2dz)(-3x1 + 4x2 - x3) + o(dz ) (left boundary,     (1)
%...    z                                         z = z0)
%...
%...                                   2
%...  xi  = (1/2dz)(xi+1 - xi-1) + o(dz ) (interior point,        (2)
%...    z                                   z ne z0, zL)
%...
%...                                          2
%...  xn  = (1/2dz)(3xn - 4xn-1 + xn-2) + o(dz ) (right boundary, (3)
%...    z                                            z = zL)
%...
%...  equations (1) to (3) apply over a grid in z with corresponding
%...  values of the function x(z) represented as
%...
%...   x1      x2       x3         xi        xn-2      xn-1    xn
%...
%...  z=z0  z=z0+dz  z=z0+2dz ... z=zi ... z=zL-2dz  z=zL-dz  z=zL
%...
%...  the origin of equations (1) to (3) is outlined below.
%...
%...  consider the following polynomial in z of arbitrary order
%...
%...                                     2             3
%...  x(z) = a0 + a1(z - z0) + a2(z - z0)  + a3(z - z0)  + ....   (4)
%...
%...  we seek the values of the coefficients a0, a1, a2, ... for a
%...  particular function x(z).  if z = z0 is substituted in equation
%...  (4), we have immediately a0 = x(z0).  next, if equation (4) is
%...  differentiated with respect to z,
%...
%...                                                   2
%...  dx(z)/dz = x (z) = a1 + 2a2(z - z*) + 3a3(z - z*)  + ...    (5)
%...              z
%...
%...  again, with z = z*, a1 = dx(z*)/dz = x (z*).  differentiation
%...                                        z
%...  of equation (5) in turn gives
%...
%...  d2x(z)/dz2 = x  (z) = 2a2 + 6a3(z - z*) + ...
%...                2z
%...
%...  and for z = z*, a2 = x  (z*)/2f (2f = 1*2, i.e., 2 factorial).
%...                        2z
%...
%...  we can continue this process of differentiation followed by the
%...  substitution z = z* to obtain the successive coefficients in
%...  equation (4), a3, a4, ...  finally, substitution of these co-
%...  efficients in equation (4) gives
%...
%...                                                 2
%...  x(z) = x(z*) + x (z*)(z - z*) + x  (z*)(z - z*)  +
%...                  z       1f       2z       2f
%...                                                              (6)
%...                                3                  4
%...                 x  (z*)(z - z*)  + x  (z*)(z - z*)  + ...
%...                  3z       3f        4z       4f
%...
%...  the correspondence between equation (6) and the well-known
%...  taylor series should be clear.  thus the expansion of a
%...  function, x(z), around a neighboring point z* in terms of x(z*)
%...  and the derivatives of x(z) at z = z* is equivalent to approxi-
%...  mating x(z) near z* by a polynomial.
%...
%...  equation (6) is the starting point for the derivation of the
%...  classical finite difference approximations of derivatives such
%...  as the three-point formulas of equations (1), (2) and (3).  we
%...  will now consider the derivation of these three-point formulas
%...  in a standard format which can then be extended to higher
%...  multi-point formulas in other functions, e.g., five-point
%...  formulas in function five_point_centered_uniform_D1.
%...
%...  three-point formulas
%...
%...     (1)  left end, point i = 1
%...
%...  if equation (6) is written around the point z = z0 for z = z0 +
%...  dz and z = z0 + 2dz, for which the corresponding values of x(z)
%...  are x1, x2 and x3 (x1 and x2 are separated with respect to z by
%...  distance dz as are x2 and x3, i.e., we assume a uniform grid
%...  spacing, dz, for independent variable z)
%...
%...                                2            3
%...  x2 = x1 + x1 ( dz) + x1  ( dz)  + x1  ( dz)  + ...          (7)
%...              z  1f      2z  2f       3z  3f
%...
%...                                2            3
%...  x3 = x1 + x1 (2dz) + x1  (2dz)  + x1  (2dz)  + ...          (8)
%...              z  1f      2z  2f       3z  3f
%...
%...  we can now take a linear combination of equations (7) and (8)
%...  by first multiplying equation (7) by a constant, a, and equa-
%...  tion (8) by constant b
%...
%...                                  2           3
%...  a(x2 = x1 + x1 ( dz) + x1  ( dz) + x1  ( dz) + ...)         (9)
%...                z  1f      2z  2f      3z  3f
%...
%...                                  2           3
%...  b(x3 = x1 + x1 (2dz) + x1  (2dz) + x1  (2dz) + ...)        (10)
%...                z  1f      2z  2f      3z  3f
%...
%...  constants a and b are then selected so that the coefficients of
%...  the x1  terms sum to one (since we are interested in obtaining
%...        z
%...  a finite difference approximation for this first derivative).
%...  also, we select a and b so that the coefficients of the x1
%...                                                            2z
%...  terms sum to zero in order to drop out the contribution of this
%...  second derivative (the basic idea is to drop out as many of the
%...  derivatives as possible in the taylor series beyond the deri-
%...  vative of interest, in this case x1 , in order to produce a
%...                                     z
%...  finite difference approximation for the derivative of maximum
%...  accuracy).  in this case we have only two constants, a and b,
%...  to select so we can drop out only the second derivative, x1  ,
%...                                                             2z
%...  in the taylor series (in addition to retaining the first deri-
%...  vative).  this procedure leads to two linear algebraic equa-
%...  tions in the two constants
%...
%...  a + 2b = 1
%...
%...  a + 4b = 0
%...
%...  solution of these equations for a and b gives
%...
%...  a = 2, b = -1/2
%...
%...  solution of equations (9) and (10) for x1  with these values of
%...  a and b gives equation (1)               z
%...
%...                                      2
%...  x1 = (1/2dz)(-3x1 + 4x2 - x3) + o(dz )                      (1)
%...    z
%...               2
%...  the term o(dz ) indicates a principal error term due to trunca-
%...                                                2
%...  tion of the taylor series which is of order dz .  this term in
%...                    2
%...  fact equals x1  dz /3f, which is easily obtained in deriving
%...                3z
%...  equation (1).
%...
%...  this same basic procedure can now be applied to the derivation
%...  of equations (2) and (3).
%...
%...     (2)  interior point i
%...
%...                                    2           3
%...  a(xi-1 = xi + xi (-dz) + xi  (-dz) + xi  (-dz) + ...)
%...                  z  1f      2z  2f      3z  3f
%...
%...                                    2           3
%...  b(xi+1 = xi + xi ( dz) + xi  ( dz) + xi  ( dz) + ...)
%...                  z  1f      2z  2f      3z  3f
%...
%...  -a + b = 1
%...
%...   a + b = 0
%...
%...  a = 1/2, b = -1/2
%...                                   2
%...  xi  = (1/2dz)(xi+1 - xi-1) + o(dz )                         (2)
%...    z
%...
%...     (3)  right end, point i = n
%...
%...                                      2            3
%...  a(xn-2 = xn + xn (-2dz) + xn  (-2dz) + xn  (-2dz) + ...)
%...                  z   1f      2z   2f      3z   3f
%...
%...                                      2            3
%...  b(xn-1 = xn + xn ( -dz) + xn  ( -dz) + xn  ( -dz) + ...)
%...                  z   1f      2z   2f      3z   3f
%...
%...  -2a - b = 1
%...
%...   4a + b = 0
%...
%...   a = -2, b = 1/2
%...                                          2
%...  xn  = (1/2dz)(3xn - 4xn-1 + xn-2) + o(dz )                  (3)
%...    z
%...
%...  the weighting coefficients for equations (1), (2) and (3) can
%...  be summarized as
%...
%...          -3   4  -1
%...
%...     1/2  -1   0   1
%...
%...           1  -4   3
%...
%...  which are the coefficients reported by bickley for n = 2, m =
%...  1, p = 0, 1, 2 (bickley, w. g., formulae for numerical differ-
%...  entiation, math. gaz., vol. 25, 1941).
%...
%...  equations (1), (2) and (3) can now be programmed to generate
%...  the derivative x (z) of x(z).
%...                  z
%...
%...  compute the spatial increment
      dz = (zL-z0)/(n-1);
      r2fdz = 1/(2*dz);
%...
%...  sparse discretization matrix      
%...
%...  interior points      
      e = ones(n,1);
      D = spdiags([-e e], [-1 1], n, n);
%...
%...  boundary points      
      D(1,1:3) = [-3 +4 -1];
%...      
      D(n,(n-2):n) = [+1 -4 +3];
%...      
      D = r2fdz*D;
