%... The MatMol Group (2016)
     function xt = chemical_reactor_1D_pdes(t,x)
%...
%... set global variables
     global D alpha rho Cp dH E k0 R;
     global cain Tin Tw v h
     global z0 zL z nz D1z D2z rL;
%...
%... transfer dependent variables
     ca=x(1:2:2*nz-1,1);
     T=x(2:2:2*nz,1);
%...
%...   Temporal derivatives
%...
%...     Entering conditions (z = 0)
         ca(1)=cain;
         cat(1)=0;
         T(1)=Tin;
         Tt(1)=0;
%...
%...     spatial derivatives
%...
         caz=D1z*ca;
         cazz= D2z*ca;
         Tz=D1z*T;
         Tzz=D2z*T;
%...
%...      Rest of reactor
          kr(2:nz,1) = k0*exp(-E./(R*T(2:nz)));
          ca_t(2:nz,1) = -v*caz(2:nz) + D*cazz(2:nz) - kr(2:nz).*ca(2:nz).^2;
          T_t(2:nz,1) = -v*Tz(2:nz) + alpha*Tzz(2:nz) + (dH/(rho*Cp))*kr(2:nz).*ca(2:nz).^2 + (2.0*h/rL)*(Tw-T(2:nz))/(rho*Cp);
%...
%... transfer temporal derivatives
     xt(1:2:2*nz-1,1)=ca_t; 
     xt(2:2:2*nz,1)=T_t;
