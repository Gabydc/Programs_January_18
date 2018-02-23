      function [D]=eleven_point_centered_uni_D1(z0,zL,n)
%...
%...  The MatMol Group (2009)
%...
%...  function eleven_point_centered_uni_D1 returns the differentiation matrix 
%...  for computing the first derivative, xz, of a variable x over the spatial domain
%...  z0 < z < zL from classical eleven-point, tenth-order finite difference 
%...  approximations (this function replaces dss0010)
%...
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
%...  the mathematical details of the finite difference approxima-
%...  tions can be summarized by the following bickley matrix
%...  for n = 10, m = 1, p = 0 to 10 (bickley, w. g., formulae for
%...  numerical differentiation, math. gaz., vol. 25, 1941)
%...
%...      -10626840            36288000           -81648000
%...      145152000          -190512000           182891520
%...
%...        -362880            -6636960            16329600
%...      -21772800            25401600           -22861440
%...
%...          40320             -806400            -4419360
%...        9676800            -8467200             6773760
%...
%...         -10080              151200            -1360800
%...       -2756160             6350400            -3810240
%...
%...           4320              -57600              388800
%...       -2073600            -1330560             4354560
%...
%...          -2880               36000             -216000
%...         864000            -3024000                   0
%...
%...           2880              -34560              194400
%...        -691200             1814400            -4354560
%...
%...          -4320               50400             -272160
%...         907200            -2116800             3810240
%...
%...          10080             -115200              604800
%...       -1935360             4233600            -6773760
%...
%...         -40320              453600            -2332800
%...        7257600           -15240960            22861440
%...
%...         362880            -4032000            20412000
%...      -62208000           127008000          -182891520
%...
%...  each entry in this table should be multiplied by 1/10f to
%...  obtain the final weighting coefficients.  from this bickley
%...  matrix, the finite difference approximation of the first
%...  derivative can be programmed for each of the grid points 1, 2,
%...  3, 4, 5,..., i,..., n-4, n-3, n-2, n-1 and n (taking into
%...  account the symmetry properties of the matrix).
%...
%...
%...  compute the spatial increment
      dz=(zL-z0)/(n-1);
      r10fdz=1./(3628800.*dz);
%...
%...  sparse discretization matrix      
%...
%...  interior points      
      e = ones(n,1);
      D = spdiags([-2880*e +36000*e -216000*e +864000*e -3024000*e 0*e +3024000*e -864000*e +216000*e -36000*e +2880*e], [-5:5], n, n);
%...
%...  boundary points      
      D([1 2 3 4 5],1:11) = [-10628640 +36288000 -81648000 +145152000 -190512000 +182891520 -127008000 +62208000 -20412000 +4032000 -362880;
                             -362880 -6636960 +16329600 -21772800 +25401600 -22861440 +15240960 -7257600 +2332800 -453600 +40320;
                             +40320 -806400 -4419360 +9676800 -8467200 +6773760 -4233600 +1935360 -604800 +115200 -10080;
                             -10080 +151200 -1360800 -2756160 +6350400 -3810240 +2116800 -907200 +272160 -50400 +4320;
                             +4320 -57600 +388800 -2073600 -1330560 +4354560 -1814400 +691200 -194400 +34560 -2880];
%...                         
      D([n-4 n-3 n-2 n-1 n],(n-10):n) = [+2880 -34560 +194400 -691200 +1814400 -4354560 +1330560 +2073600 -388800 +57600 -4320;
                                         -4320 +50400 -272160 +907200 -2116800 +3810240 -6350400 +2756160 +1360800 -151200 +10080;
                                         +10080 -115200 +604800 -1935360 +4233600 -6773760 +8467200 -9676800 +4419360 +806400 -40320;
                                         -40320 +453600 -2332800 +7257600 -15240960 +22861440 -25401600 +21772800 -16329600 +6636960 +362880; 
                                         +362880 -4032000 +20412000 -62208000 +127008000 -182891520 +190512000 -145152000 +81648000 -36288000 +10628640];
%...      
      D=r10fdz*D;



         

 
 