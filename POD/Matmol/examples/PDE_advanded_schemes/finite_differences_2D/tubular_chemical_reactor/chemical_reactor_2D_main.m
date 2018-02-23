%... The MatMol Group (2016)
%... Dynamic analysis of a 2-D reactor
%...
%... The following equations model a 2-D reactor
%...
%...   Material balance
%...
%...   ca  = -v*ca  + D*(ca + ca   + (1/r)ca ) - kr*ca^2                (1)
%...     t        z        zz   rr          r
%...
%...   Energy balance
%...
%...   T  = -v*T  + lambda*/(rho*Cp)*(T + T   + (1/r)*T )              (2)
%...    t       z                      zz  rr          r
%...
%...        - dH*kr*ca^2/(rho*Cp)
%...
%...   kr = k0*exp(-E/(R*T))                                       (3)
%...
%...   T  (rL,z,t) = h*(T  - T(rL,z,t))                            (4)
%...    r                w
%...
%... The variables and parameters for this model are (in cgs units)
%...
%...   Reactant concentration          ca        (eq. (1))
%...
%...   Temperature                     T         (eq. (3))
%...
%...   Reaction rate constant          kr        (eq. (3))
%...
%...   Time                            t
%...
%...   Radial position                 r
%...
%...   Axial position                  z
%...
%...   Entering concentration          cain         0.01
%...
%...   Entering temperature            Tin          305
%...
%...   Wall temperature                Tw           355
%...
%...   Reactor radius                  rL           2
%...
%...   Reactor length                  zL           100
%...
%...   Maximum linear velocity         vmax         1
%...
%...   Mass diffusivity                D            0.1
%...
%...   Thermal diffusivity  alpha=lambda/(rho*Cp)   0.1
%...
%...   Liquid density                  rho          1.0
%...
%...   Liquid specific heat            Cp           0.5
%...
%...   Heat of reaction                -dH          10000
%...                                                       
%...   Specific rate constant          k0           2.0e+09
%...
%...   Activation energy               E            15000
%...
%...   Gas constant                    R            1.987
%...
%...   heat transfert with jacket      h            0.05
%...
%... Two cases can be considered:
%...
%...    (1)  Constant velocity profile, v(r) = vavg = vmax/2
%...
%...    (2)  Parabolic velocity profile with v(r) = 
%...         vmax*(1 - (r/rL)^2)
%...
%... The purpose of these two runs is to determine the effect
%... of a velcoity profile on the yield of the reactor.
%...
%... The solution to this system, ca(r,z,t) (from eq. (1)) and
%... T(r,z,t) (from eq. (2)) is computed by a fixed step Euler
%... integration in t.  Also, the choice of an Euler step is
%... sensitive because of the nonlinearity of eq. (3).
%...
     close all
     clear all
%...
%... start a stopwatch timer
     tic
%...
%... set global variables
     global D alpha rho Cp dH E k0 R;
     global cain Tin Tw v h gamma
     global r0 rL r drs z0 zL z nr nz D1r D1z D2r D2z;
%...
%... model parameters (defined in the comments above)
     Tw=355.0;
     h=0.05;
     gamma=100;
     vmax=1.0;
     D=0.1;
     alpha=0.1;
     rho=1.0;
     Cp=0.5;
     dH=10000.0;
     E=15000.0;
     k0=2.0e+09;
     R=1.987;
%...
%... inlet concentration and temperature
     cain=0.01;
     Tin=305.0;
%...
%... grid in axial direction
     z0=0.0;
     zL=100.0;
     nz=101;
     dz=(zL-z0)/(nz-1);
     z=[z0:dz:zL]';
%...
%... grid in radial direction
     r0=0.0;
     rL=2.0;
     nr=11;
     dr=(rL-r0)/(nr-1);
     r=[r0:dr:rL]';
     drs=dr^2;
%...
%... velocity profile
     v = (vmax/2)*ones(nr,1);
%      v = vmax*(1- (r/rL).^2);
%...
%... initial conditions
     ca=cain*ones(nr,nz);
     T=Tin*ones(nr,nz);
%...
%... create 1D arrays
     ca=reshape(ca,nr*nz,1);
     T=reshape(T,nr*nz,1);
     x(1:2:2*nr*nz-1,1)=ca; 
     x(2:2:2*nr*nz,1)=T;
%...
%... differentiation matrices in r
     D1r=three_point_centered_D1(r);
     D2r=three_point_centered_D2(r);
%...
%... differentiation matrix in z
     D1z=two_point_upwind_D1(z,vmax);
     D2z=three_point_centered_D2(z);
%...  
%... call to ODE solver 
%...
     t0=0.0;
     tf=300.0;
     Dt=0.1;
     Dtplot=20.0;
     options=[];
%...     
     [tout, yout] = ode15s(@chemical_reactor_2D_pdes,[t0:Dtplot:tf],x,options);
%     [tout, yout] = euler_solver(@reactor_2D_pdes,t0,tf,x,Dt,Dtplot);
%...     
     for k=1:length(tout),
         ca_out=reshape(yout(k,1:2:2*nr*nz-1),nr,nz);
         figure(1);
         surf(z,r,ca_out,'EdgeColor','none');
         xlabel('z');
         ylabel('r');
         zlabel('ca(r,z,t)');
         title('Concentration')
         axis([0 100 0 2 0 0.01]);
         pause(1)
     end

     for k=1:length(tout),
         T_out=reshape(yout(k,2:2:2*nr*nz),nr,nz);
         figure(2);
         surf(z,r,T_out,'EdgeColor','none');
         xlabel('z');
         ylabel('r');
         zlabel('T(r,z,t)');
         title('Temperature')
         axis([0 100 0 2 300 550]);
         pause(1)
     end
%...
%... read the stopwatch timer
     tcpu=toc;
