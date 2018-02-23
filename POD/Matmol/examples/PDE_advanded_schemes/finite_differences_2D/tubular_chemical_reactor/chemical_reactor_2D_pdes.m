%... The MatMol Group (2016)
     function xt = chemical_reactor_2D_pdes(t,x)
%...
%... set global variables
     global D alpha rho Cp dH E k0 R;
     global cain Tin Tw v h gamma
     global r0 rL r drs z0 zL z nr nz D1r D1z D2r D2z;
     
%...
%... create 2D arrays from x
     ca=reshape(x(1:2:2*nr*nz-1,1),nr,nz);
     T=reshape(x(2:2:2*nr*nz,1),nr,nz);
%...
%...   Temporal derivatives
%...
%...     Cover the nr x nz grid
%...
%...       Entering conditions (z = 0)
           for i=1:nr
             ca(i,1)=cain;
             cat(i,1)=0;
             T(i,1)=Tin;
             Tt(i,1)=0;
           end
%...
%...     spatial derivatives
%...
%...     car  (by three-point centered approximations)
         nd=1;
         car=diff_2D(nr,nz,nd,D1r,ca);
%...
%...     carr  (by three-point centered approximations)
         nd=1;
         carr=diff_2D(nr,nz,nd,D2r,ca);
%...
%...     caz  (by two-point upwind approximations)
         nd=2;
         caz=diff_2D(nr,nz,nd,D1z,ca);
%...
%...     cazz  (by three-point centered approximations)
         nd=2;
         cazz=diff_2D(nr,nz,nd,D2z,ca);
%...
%...     Tr  (by three-point centered approximations)
         nd=1;
         Tr=diff_2D(nr,nz,nd,D1r,T);
%...
%...     Trr  (by three-point centered approximations)
         nd=1;
         Trr=diff_2D(nr,nz,nd,D2r,T);
%...
%...     Tz  (by two-point upwind approximations)
         nd=2;
         Tz=diff_2D(nr,nz,nd,D1z,T);
%...
%...     Tzz  (by three-point centered approximations)
         nd=2;
         Tzz=diff_2D(nr,nz,nd,D2z,T);
%...
%...       Rest of reactor
           for j=2:nz
%...
%...         Centerline (r = 0)
             kr=k0*exp(-E/(R*T(1,j)));
             ca_t(1,j)=-v(1)*caz(1,j)+4.0*D*(ca(2,j)-ca(1,j))/drs-kr*ca(1,j)^2;
             T_t(1,j)=-v(1)*Tz(1,j)+4.0*alpha*(T(2,j)-T(1,j))/drs+dH*kr/(rho*Cp)*ca(1,j)^2;
%...   
%...         Wall (r = rL)
             kr=k0*exp(-E/(R*Tw));
             ca_t(nr,j)=-v(nr)*caz(nr,j)+2.0*D*(ca(nr-1,j)-ca(nr,j))/drs-kr*ca(nr,j)^2;
             T_t(nr,j)=gamma*(h*(Tw-T(nr,j))-Tr(nr,j));
%...
%...         Interior, r ~= 0 and r ~= rL
             for i=2:nr-1
%...
             kr=k0*exp(-E/(R*T(i,j)));
             ca_t(i,j)=-v(i)*caz(i,j)+D*(carr(i,j)+(1/r(i))*car(i))-kr*ca(i,j)^2;
             T_t(i,j)=-v(i)*Tz(i,j)+alpha*(Trr(i,j)+(1/r(i))*Tr(i,j))+dH*kr/(rho*Cp)*ca(i,j)^2;
%...
%...         next r
             end
%...
%...       next z
           end

%...
%... create a 1D array from cat and Tt
     ca_t=reshape(ca_t,nr*nz,1);
     T_t=reshape(T_t,nr*nz,1);
     xt(1:2:2*nr*nz-1,1)=ca_t; 
     xt(2:2:2*nr*nz,1)=T_t;
