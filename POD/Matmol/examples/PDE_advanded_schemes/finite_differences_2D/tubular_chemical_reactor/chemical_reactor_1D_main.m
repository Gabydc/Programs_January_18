%... The MatMol Group (2016)
%... Dynamic analysis of a 1-D reactor
%...
%... The following equations model a 1-D reactor
%...
%...   Material balance
%...
%...   ca  = -v*ca  + D*ca   - kr*ca^2                                        (1)
%...     t        z       zz
%...
%...   Energy balance
%...
%...   T  = -v*T  + lambda/(rho*Cp)*T   + (-dH)*kr*ca^2/(rho*Cp) + (2*h/rL) (T  - T)/(rho*Cp) (2)
%...    t       z                    zz                                    w
%...
%...   kr = k0*exp(-E/(R*T))                                          (3)
%...
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
%...   Axial position                  z
%...
%...   Entering concentration          cain         0.01
%...
%...   Entering temperature            Tin          305
%...
%...   Wall temperature                Tw           355
%...
%...   Reactor length                  zL           100
%...
%...   Maximum linear velocity         vmax         1
%...
%...   Mass diffusivity                D            0.1
%...
%...   Thermal diffusivity alpha = lambda/(rho*Cp)   0.1
%...
%...   Liquid density                  rho          1.0
%...
%...   Liquid specific heat            Cp           0.5
%...
%...   Heat of reaction                dH          -10000
%...                                                       
%...   Specific rate constant          k0           2.0e+09
%...
%...   Activation energy               E            15000
%...
%...   Gas constant                    R            1.987
%...
%...   heat transfert with jacket      h            0.05
%...
%...
%... The purpose of these two runs is to determine the effect
%... of a velcity profile on the yield of the reactor.
%...
%... The solution to this system, ca(z,t) (from eq. (1)) and
%... T(z,t) (from eq. (2)) is computed by a fixed step Euler
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
     global cain Tin Tw v h
     global z0 zL z nz D1z D2z rL;
%...
%... model parameters (defined in the comments above)
     Tw=355.0;
     h=0.05;
     rL=2;
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
%... velocity
     v = vmax/2;
%...
%... initial conditions
     ca=cain*ones(nz,1);
     T=Tin*ones(nz,1);
%...
%... transfer dependent variables
     x(1:2:2*nz-1)=ca; 
     x(2:2:2*nz)=T;
%...
%... differentiation matrix in z (convective and diffusive terms)
     D1z=two_point_upwind_D1(z,v);
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
     [tout, yout] = ode15s(@chemical_reactor_1D_pdes,[t0:Dtplot:tf],x,options);
%     [tout, yout] = euler_solver(@reactor_1D_pdes,t0,tf,x,Dt,Dtplot);
%...     
     for k=1:length(tout),
         ca_out=yout(k,1:2:2*nz-1);
         figure(1);
         hold on
         plot(z,ca_out);
         xlabel('z');
         ylabel('ca(z,t)');
         title('Concentration')
         axis([0 100 0 0.01]);
     end
     for k=1:length(tout),
         T_out=yout(k,2:2:2*nz);
         figure(2);
         hold on
         plot(z,T_out);
         xlabel('z');
         ylabel('T(z,t)');
         title('Temperature')
         axis([0 100 300 450]);
     end     
%...
%... read the stopwatch timer
     tcpu=toc;
