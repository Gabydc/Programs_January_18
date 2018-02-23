%... The MatMol Group (2016)
     function xt = bioreactor_pdes(t,x)
%...
%... set global variables
     global v nu1 eps mumax Kc kd;
     global Sin gamma
     global z0 zL z nz D1z;
%...
%... transfer dependent variables
     S = x(1:2:2*nz-1,1);
     X = x(2:2:2*nz,1);
%...
%...   Temporal derivatives
%...
%...     Entering conditions (z = 0)
         S_t(1,1) = gamma*(Sin-S(1));
%...
%...     kinetics
         phi1 = mumax*S.*X./(Kc*X+S);
         phi2 = kd*X;
%...
%...     spatial derivatives
%...
%...     Sz  (by two-point upwind approximations)
         Sz = D1z*S;
%...
%...      Rest of bioreactor
          S_t(2:nz,1) = -v*Sz(2:nz) - nu1*(1-eps)*phi1(2:nz)/eps;
          X_t = phi1 - phi2;
%...
%... transfer temporal derivatives
     xt(1:2:2*nz-1,1) = S_t; 
     xt(2:2:2*nz,1) = X_t;
