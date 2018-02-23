      function [D] = five_point_centered_uni_D1(z0,zL,n)
%...
%...  The MatMol Group (2009)
%...
%...  function five_point_centered_uni_D1 returns the differentiation matrix 
%...  for computing the first derivative, xz, of a variable u over the spatial domain
%...  z0 < z < zL from classical five-point, fourth-order finite difference 
%...  approximations (this function replaces dss004)
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
%...  five-point formulas
%...
%...     (1)  left end, point i = 1
%...
%...                                   2            3            4
%...  a(x2 = x1 + x1  ( dz) + x1  ( dz)  + x1  ( dz)  + x1  ( dz)
%...                z   1f      2z  2f       3z  3f       4z  4f
%...
%...                      5            6            7
%...           + x1  ( dz)  + x1  ( dz)  + x1  ( dz)  + ...)
%...               5z  5f       6z  6f       7z  7f
%...
%...                                   2            3            4
%...  b(x3 = x1 + x1  (2dz) + x1  (2dz)  + x1  (2dz)  + x1  (2dz)
%...                z   1f      2z  2f       3z  3f       4z  4f
%...
%...                      5            6            7
%...           + x1  (2dz)  + x1  (2dz)  + x1  (2dz)  + ...)
%...               5z  5f       6z  6f       7z  7f
%...
%...                                   2            3            4
%...  c(x4 = x1 + x1  (3dz) + x1  (3dz)  + x1  (3dz)  + x1  (3dz)
%...                z   1f      2z  2f       3z  3f       4z  4f
%...
%...                      5            6            7
%...           + x1  (3dz)  + x1  (3dz)  + x1  (3dz)  + ...)
%...               5z  5f       6z  6f       7z  7f
%...
%...                                   2            3            4
%...  d(x5 = x1 + x1  (4dz) + x1  (4dz)  + x1  (4dz)  + x1  (4dz)
%...                z   1f      2z  2f       3z  3f       4z  4f
%...
%...                      5            6            7
%...           + x1  (4dz)  + x1  (4dz)  + x1  (4dz)  + ...)
%...               5z  5f       6z  6f       7z  7f
%...
%...  constants a, b, c and d are selected so that the coefficients
%...  of the x1  terms sum to one and the coefficients of the x1  ,
%...           z                                                2z
%...  x1   and x1   terms sum to zero
%...    3z       4z
%...
%...  a +   2b +   3c +   4d = 1
%...
%...  a +   4b +   9c +  16d = 0
%...
%...  a +   8b +  27c +  64d = 0
%...
%...  a +  16b +  81c + 256d = 0
%...
%...  simultaneous solution for a, b, c and d followed by the solu-
%...  tion of the preceding taylor series, truncated after the x
%...                                                            4z
%...  terms, for x1  gives the following five-point approximation
%...               z
%...                                                         4
%...  x1  = (1/12dz)(-25x1 + 48x2 - 36x3 + 16x4 - 3x5) + o(dz )   (1)
%...    z
%...
%...     (2)  interior point, i = 2
%...
%...                                   2            3            4
%...  a(x1 = x2 + x2  (-dz) + x2  (-dz)  + x2  (-dz)  + x2  (-dz)
%...                z   1f      2z  2f       3z  3f       4z  4f
%...
%...                      5            6            7
%...           + x2  (-dz)  + x2  (-dz)  + x2  (-dz)  + ...)
%...               5z  5f       6z  6f       7z  7f
%...
%...                                   2            3            4
%...  b(x3 = x2 + x2  ( dz) + x2  ( dz)  + x2  ( dz)  + x2  ( dz)
%...                z   1f      2z  2f       3z  3f       4z  4f
%...
%...                      5            6            7
%...           + x2  ( dz)  + x2  ( dz)  + x2  ( dz)  + ...)
%...               5z  5f       6z  6f       7z  7f
%...
%...                                   2            3            4
%...  c(x4 = x2 + x2  (2dz) + x2  (2dz)  + x2  (2dz)  + x2  (2dz)
%...                z   1f      2z  2f       3z  3f       4z  4f
%...
%...                      5            6            7
%...           + x2  (2dz)  + x2  (2dz)  + x2  (2dz)  + ...)
%...               5z  5f       6z  6f       7z  7f
%...
%...                                   2            3            4
%...  d(x5 = x2 + x2  (3dz) + x2  (3dz)  + x2  (3dz)  + x2  (3dz)
%...                z   1f      2z  2f       3z  3f       4z  4f
%...
%...                      5            6            7
%...           + x2  (3dz)  + x2  (3dz)  + x2  (3dz)  + ...)
%...               5z  5f       6z  6f       7z  7f
%...
%...  -a +   b +  2c +  3d = 1
%...
%...   a +   b +  4c +  9d = 0
%...
%...  -a +   b +  8c + 27d = 0
%...
%...   a +   b + 16c + 81d = 0
%...
%...  simultaneous solution for a, b, c and d followed by the solu-
%...  tion of the preceding taylor series, truncated after the x
%...                                                            4z
%...  terms, for x1  gives the following five-point approximation
%...               z
%...                                                        4
%...  x2  = (1/12dz)(-3x1 - 10x2 + 18x3 -  6x4 +  x5) + o(dz )    (2)
%...    z
%...
%...     (3)  interior point i, i ne 2, n-1
%...
%...                                        2             3
%...  a(xi-2 = xi + xi  (-2dz)  + xi  (-2dz)  + xi  (-2dz)
%...                  z    1f       2z   2f       3z   3f
%...
%...                          4             5             6
%...              + xi  (-2dz)  + xi  (-2dz)  + xi  (-2dz)  + ...)
%...                  4z   4f       5z   5f       6z   6f
%...
%...                                        2             3
%...  b(xi-1 = xi + xi  ( -dz)  + xi  ( -dz)  + xi  ( -dz)
%...                  z    1f       2z   2f       3z   3f
%...
%...                          4             5             6
%...              + xi  ( -dz)  + xi  ( -dz)  + xi  ( -dz)  + ...)
%...                  4z   4f       5z   5f       6z   6f
%...
%...                                        2             3
%...  c(xi+1 = xi + xi  (  dz)  + xi  (  dz)  + xi  (  dz)
%...                  z    1f       2z   2f       3z   3f
%...
%...                          4             5             6
%...              + xi  (  dz)  + xi  (  dz)  + xi  (  dz)  + ...)
%...                  4z   4f       5z   5f       6z   6f
%...
%...                                        2             3
%...  d(xi+2 = xi + xi  ( 2dz)  + xi  ( 2dz)  + xi  ( 2dz)
%...                  z    1f       2z   2f       3z   3f
%...
%...                          4             5             6
%...              + xi  ( 2dz)  + xi  ( 2dz)  + xi  ( 2dz)  + ...)
%...                  4z   4f       5z   5f       6z   6f
%...
%...   -2a -   b +   c +  2d = 1
%...
%...    4a +   b +   c +  4d = 0
%...
%...   -8a -   b +   c +  8d = 0
%...
%...   16a +   b +   c + 16d = 0
%...
%...  simultaneous solution for a, b, c and d followed by the solu-
%...  tion of the preceding taylor series, truncated after the x
%...                                                            4z
%...  terms, for x1  gives the following five-point approximation
%...               z
%...                                                          4
%...  xi  = (1/12dz)(xi-2 - 8xi-1 + 0xi + 8xi+1 - xi+2) + o(dz )  (3)
%...    z
%...
%...     (4)  interior point, i = n-1
%...
%...                                              2               3
%...  a(xn-4 = xn-1 + xn-1  (-3dz)  + xn-1  (-3dz)  + xn-1  (-3dz)
%...                      z    1f         2z   2f         3z   3f
%...
%...                       4               5               6
%...         + xn-1  (-3dz)  + xn-1  (-3dz)  + xn-1  (-3dz)  + ...
%...               4z   4f         5z   5f         6z   6f
%...
%...                                              2               3
%...  b(xn-3 = xn-1 + xn-1  (-2dz)  + xn-1  (-2dz)  + xn-1  (-2dz)
%...                      z    1f         2z   2f         3z   3f
%...
%...                       4               5               6
%...         + xn-1  (-2dz)  + xn-1  (-2dz)  + xn-1  (-2dz)  + ...
%...               4z   4f         5z   5f         6z   6f
%...
%...                                              2               3
%...  c(xn-2 = xn-1 + xn-1  ( -dz)  + xn-1  ( -dz)  + xn-1  ( -dz)
%...                      z    1f         2z   2f         3z   3f
%...
%...                       4               5               6
%...         + xn-1  ( -dz)  + xn-1  ( -dz)  + xn-1  ( -dz)  + ...
%...               4z   4f         5z   5f         6z   6f
%...
%...                                              2               3
%...  d(xn   = xn-1 + xn-1  (  dz)  + xn-1  (  dz)  + xn-1  (  dz)
%...                      z    1f         2z   2f         3z   3f
%...
%...                       4               5               6
%...         + xn-1  (  dz)  + xn-1  (  dz)  + xn-1  (  dz)  + ...
%...               4z   4f         5z   5f         6z   6f
%...
%...  -3a -  2b -   c +   d = 1
%...
%...   9a +  4b +   c +   d = 0
%...
%... -27a -  8b -   c +   d = 0
%...
%...  81a + 16b +   c +   d = 0
%...
%...  simultaneous solution for a, b, c and d followed by the solu-
%...  tion of the preceding taylor series, truncated after the x
%...                                                            4z
%...  terms, for x1  gives the following five-point approximation
%...               z
%...                                                                4
%...  xn-1  = (1/12dz)(-xn-4 + 6xn-3 - 18xn-2 + 10xn-1 + 3xn) + o(dz )
%...      z
%...                                                              (4)
%...
%...    (5)  right end, point i = n
%...
%...                                       2             3
%...  a(xn-4 = xn + xn (-4dz)  + xn  (-4dz)  + xn  (-4dz)
%...                  z   1f       2z   2f       3z   3f
%...
%...                         4             5             6
%...             + xn  (-4dz)  + xn  (-4dz)  + xn  (-4dz)  + ...)
%...                 4z   4f       5z   5f       6z   6f
%...
%...                                       2             3
%...  b(xn-3 = xn + xn (-3dz)  + xn  (-3dz)  + xn  (-3dz)
%...                  z   1f       2z   2f       3z   3f
%...
%...                         4             5             6
%...             + xn  (-3dz)  + xn  (-3dz)  + xn  (-3dz)  + ...)
%...                 4z   4f       5z   5f       6z   6f
%...
%...                                       2             3
%...  c(xn-2 = xn + xn (-2dz)  + xn  (-2dz)  + xn  (-2dz)
%...                  z   1f       2z   2f       3z   3f
%...
%...                         4             5             6
%...             + xn  (-2dz)  + xn  (-2dz)  + xn  (-2dz)  + ...)
%...                 4z   4f       5z   5f       6z   6f
%...
%...                                       2             3
%...  d(xn-1 = xn + xn ( -dz)  + xn  ( -dz)  + xn  ( -dz)
%...                  z   1f       2z   2f       3z   3f
%...
%...                         4             5             6
%...             + xn  ( -dz)  + xn  ( -dz)  + xn  ( -dz)  + ...)
%...                 4z   4f       5z   5f       6z   6f
%...
%...   -4a -  3b -  2c -   d = 1
%...
%...   16a +  9b +  4c +   d = 0
%...
%...  -64a - 27b -  8c -   d = 0
%...
%...  256a + 81b + 16c +   d = 0
%...
%...  simultaneous solution for a, b, c and d followed by the solu-
%...  tion of the preceding taylor series, truncated after the x
%...                                                            4z
%...  terms, for x1  gives the following five-point approximation
%...               z
%...                                                                4
%...  xn  = (1/12dz)(3xn-4 - 16xn-3 + 36xn-2 - 48xn-1 + 25xn) + o(dz )
%...    z
%...                                                              (5)
%...
%...  the weighting coefficients for equations (1) to (5) can be
%...  summarized as
%...
%...             -25   48  -36   16   -3
%...
%...              -3  -10   18   -6    1
%...
%...       1/12    1   -8    0    8   -1
%...
%...              -1    6  -18   10    3
%...
%...               3  -16   36  -48   25
%...
%...  which are the coefficients reported by bickley for n = 4, m =
%...  1, p = 0, 1, 2, 3, 4 (bickley, w. g., formulae for numerical
%...  differentiation, math. gaz., vol. 25, 1941.  note - the bickley
%...  coefficients have been divided by a common factor of two).
%...
%...  equations (1) to (5) can now be programmed to generate the
%...  derivative x (z) of function x(z).
%...              z
%...
%...  compute the spatial increment
      dz=(zL-z0)/(n-1);
      r4fdz=1/(12*dz);
%...
%...  sparse discretization matrix      
%...
%...  interior points      
      e = ones(n,1);
      D = spdiags([+1*e -8*e 0*e +8*e -1*e], [-2:2], n, n);
%...
%...  boundary points      
      D([1 2],1:5) = [-25 +48 -36 +16 -3; -3 -10 +18 -6 1];
%...      
      D([n-1 n],(n-4):n) = [-1 +6 -18 +10 +3; +3 -16 +36 -48 +25];
%...      
      D=r4fdz*D;
      
      
