      function [D] = nine_point_centered_unif_D1(z0,zL,n)
%...
%...  The MatMol Group (2009)
%...
%...  function nine_point_centered_unif_D1 returns the differentiation matrix 
%...  for computing the first derivative, xz, of a variable x over the spatial domain
%...  z0 < z < zL from classical nine-point, eighth-order finite difference 
%...  approximations (this function replaces dss008)
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
%...  nine-point formulas
%...
%...     (1)  left end, point i = 1
%...
%...                                  2            3            4
%...  a(x2 = x1 + x1 ( dz) + x1  ( dz)  + x1  ( dz)  + x1  ( dz)
%...                z  1f      2z  2f       3z  3f       4z  4f
%...
%...                      5            6            7
%...           + x1  ( dz)  + x1  ( dz)  + x1  ( dz)  + ...)
%...               5z  5f       6z  6f       7z  7f
%...
%...                                  2            3            4
%...  b(x3 = x1 + x1 (2dz) + x1  (2dz)  + x1  (2dz)  + x1  (2dz)
%...                z  1f      2z  2f       3z  3f       4z  4f
%...
%...                      5            6            7
%...           + x1  (2dz)  + x1  (2dz)  + x1  (2dz)  + ...)
%...               5z  5f       6z  6f       7z  7f
%...
%...                                  2            3            4
%...  c(x4 = x1 + x1 (3dz) + x1  (3dz)  + x1  (3dz)  + x1  (3dz)
%...                z  1f      2z  2f       3z  3f       4z  4f
%...
%...                      5            6            7
%...           + x1  (3dz)  + x1  (3dz)  + x1  (3dz)  + ...)
%...               5z  5f       6z  6f       7z  7f
%...
%...                                  2            3            4
%...  d(x5 = x1 + x1 (4dz) + x1  (4dz)  + x1  (4dz)  + x1  (4dz)
%...                z  1f      2z  2f       3z  3f       4z  4f
%...
%...                      5            6            7
%...           + x1  (4dz)  + x1  (4dz)  + x1  (4dz)  + ...)
%...               5z  5f       6z  6f       7z  7f
%...
%...                                  2            3            4
%...  e(x6 = x1 + x1 (5dz) + x1  (5dz)  + x1  (5dz)  + x1  (5dz)
%...                z  1f      2z  2f       3z  3f       4z  4f
%...
%...                      5            6            7
%...           + x1  (5dz)  + x1  (5dz)  + x1  (5dz)  + ...)
%...               5z  5f       6z  6f       7z  7f
%...                                  2            3            4
%...  f(x7 = x1 + x1 (6dz) + x1  (6dz)  + x1  (6dz)  + x1  (6dz)
%...                z  1f      2z  2f       3z  3f       4z  4f
%...
%...                      5            6            7
%...           + x1  (6dz)  + x1  (6dz)  + x1  (6dz)  + ...)
%...               5z  5f       6z  6f       7z  7f
%...
%...                                  2            3            4
%...  g(x8 = x1 + x1 (7dz) + x1  (7dz)  + x1  (7dz)  + x1  (7dz)
%...                z  1f      2z  2f       3z  3f       4z  4f
%...
%...                      5            6            7
%...           + x1  (7dz)  + x1  (7dz)  + x1  (7dz)  + ...)
%...               5z  5f       6z  6f       7z  7f
%...
%...                                  2            3            4
%...  h(x9 = x1 + x1 (8dz) + x1  (8dz)  + x1  (8dz)  + x1  (8dz)
%...                z  1f      2z  2f       3z  3f       4z  4f
%...
%...                      5            6            7
%...           + x1  (8dz)  + x1  (8dz)  + x1  (8dz)  + ...)
%...               5z  5f       6z  6f       7z  7f
%...
%...  constants a, b, c, d, e, f, g and h are selected so that the
%...  coefficients of the x1  terms sum to one and the coefficients of
%...                        z
%...  the x1  , x1  , x1  , x1  , x1  , x1   and x1   sum to zero.
%...        2z    3z    4z    5z    6z,   7z       8z
%...
%...        1      1      1      1      1      1      1
%...  a +  2 b +  3 c +  4 d +  5 e +  6 f +  7 g +  8 h = 1
%...
%...        2      2      2      2      2      2      2
%...  a +  2 b +  3 c +  4 d +  5 e +  6 f +  7 g +  8 h = 0
%...
%...        3      3      3      3      3      3      3
%...  a +  2 b +  3 c +  4 d +  5 e +  6 f +  7 g +  8 h = 0
%...
%...        4      4      4      4      4      4      4
%...  a +  2 b +  3 c +  4 d +  5 e +  6 f +  7 g +  8 h = 0
%...
%...        5      5      5      5      5      5      5
%...  a +  2 b +  3 c +  4 d +  5 e +  6 f +  7 g +  8 h = 0
%...
%...        6      6      6      6      6      6      6
%...  a +  2 b +  3 c +  4 d +  5 e +  6 f +  7 g +  8 h = 0
%...
%...        7      7      7      7      7      7      7
%...  a +  2 b +  3 c +  4 d +  5 e +  6 f +  7 g +  8 h = 0
%...
%...        8      8      8      8      8      8      8
%...  a +  2 b +  3 c +  4 d +  5 e +  6 f +  7 g +  8 h = 0
%...
%...  simultaneous solution for a, b, c, d, e, f, g and h followed by
%...  the solution of the preceding taylor series, truncated after the
%...  x   terms, for x1  gives the following nine-point approximation
%...   8z              z
%...
%...  x1  = (1/8f)(-109584x1 + 322560x2 - 564480x3 + 752640x4
%...    z
%...               -705600x5 + 451584x6 - 188160x7 + 46080x8
%...                               8
%...               - 5040x9) + o(dz )                             (1)
%...
%...  the preceding analysis can be repeated to produce nine-point
%...  approximations for the first derivatives x2 , x3 , x4 , xi ,
%...                                             z    z    z    z
%...  xn-3 , xn-2 , xn-1  and xn .  the results can be summarized by
%...      z      z      z       z
%...  the following bickley matrix for n = 8, m = 1, p = 0 to 8,
%...  (bickley, w. g., formulae for numerical differentiation, math.
%...  gaz., vol. 25, 1941)
%...
%...            -109584   322560  -564480   752640  -705600
%...
%...              -5040   -64224   141120  -141120   117600
%...
%...                720   -11520   -38304    80640   -50400
%...
%...               -240     2880   -20160   -18144    50400
%...
%...     1/8f       144    -1536     8064   -32256        0
%...
%...               -144     1440    -6720    20160   -50400
%...
%...                240    -2304    10080   -26880    50400
%...
%...               -720     6720   -28224    70560  -117600
%...
%...               5040   -46080   188160  -451584   705600
%...
%...  from this bickley matrix, the finite difference approximation of
%...  the first derivative can be programmed for each of the grid points
%...  1, 2, 3, 4,..., i,..., n-3, n-2, n-1, n (taking into account
%...  the symmetry properties of the bickley matrix).
%...
%...
%...  compute the spatial increment
      dz=(zL-z0)/(n-1);
      r8fdz=1/(40320*dz);
%...
%...  sparse discretization matrix      
%...
%...  interior points      
      e = ones(n,1);
      D = spdiags([+144*e -1536*e +8064*e -32256*e 0*e +32256*e -8064*e +1536*e -144*e], [-4:4], n, n);
%...
%...  boundary points      
      D([1 2 3 4],1:9) = [-109584 +322560 -564480 +752640 -705600 +451584 -188160 +46080 -5040;
                          -5040 -64224 +141120 -141120 +117600 -70560 +28224 -6720 +720;
                          +720 -11520 -38304 +80640 -50400 +26880 -10080 +2304 -240;
                          -240 +2880 -20160 -18144 +50400 -20160 +6720 -1440 +144];
%...
      D([n-3 n-2 n-1 n],(n-8):n) = [-144 +1440 -6720 +20160 -50400 +18144 +20160 -2880 +240;
                                    +240 -2304 +10080 -26880 +50400 -80640 +38304 +11520 -720;
                                    -720 +6720 -28224 +70560 -117600 +141120 -141120 +64224 +5040; 
                                    +5040 -46080 +188160 -451584 +705600 -752640 +564480 -322560 +109584];
%...      
      D=r8fdz*D;
