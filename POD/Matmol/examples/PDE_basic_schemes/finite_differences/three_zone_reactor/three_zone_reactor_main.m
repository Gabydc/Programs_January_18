%... The MatMol Group (2016)
%...
%... Fixed bed reactor with inert portions experiencing catalyst decay
%...
%... Consider a fixed-bed reactor with three sections: a inert entrance section,
%... a catalyst section and a inert oulet section
%...
%... The reaction under consideration is benzene hydrogenation (exothermic reaction)
%... The poisoning kinetics of thiophene on Ni-kieselguhr catalyst and the deactivation
%... behavior of nonisothermal fixed-bed is investigated (Weng, Eigenberger and Butt, 1975)
%...
%... The system can be described by the following mass and energy balances:
%...
%... cBt = -v*cBz + DB*cBzz - rhoc*rB/eps
%...
%... cTt = -v*cTz + DT*cTzz - rhoc*rT/eps
%...
%... Tt = -((eps*v*rhog*cpg)/rhocp)*Tz + (leff/rhocp)*Tzz + 2*alpha*(Tw - T)/(Rr*rhocp) + (DH/rhocp)*rhoc*rB
%...
%... thetat = - rd;
%...
%... where
%...
%...  cB        benzene concentration
%...
%...  cB        thiophene concentration
%...
%...  T         temperature
%...
%...  theta     catalyst activity 
%...
%...  t         time
%...
%...  z         space
%...
%...
%... the inlet and outlet boundary conditions (BCs) are taken as Neumann BCs
%...
%... cz(z=z0,t) = (v/D)*(c(z=z0,t)-cin)
%...
%... Tz(z=z0,t) = (eps*v*rhog*cpg/leff)*(T(z=z0,t)-Tin)
%...
%... cz(z=zL,t) = 0             Tz(z=zL,t) = 0                
%...
%... where
%...
%...  [z0,zL]      spatial domain
%...
%... Reference
%...
%... Weng H.S., Eigenberger G. and Butt J.B.
%... Catalyst poisoning and fixed bed reactor dynamics
%... Chem. Engng Sci. 1975
%...
%... the following code computes a solution to this problem
%...
     close all
     clear all
%...
%... start a stopwatch timer
     tic
%...
%... set global variables
     global v DB DT rhocp cpg rhog Tw Rr DH lH2 Dp leff eps alpha rhoc;
     global ki0 K0 Q EB R_kcal R_J kd0 Ed MT;
     global Ptot cBin cTin Tin cHin xHin;
     global z01 zL1 z02 zL2 z03 zL3 n1 n2 n3 z1 z2 z3 D1_1 D2_1 D1_2 D2_2 D1_3 D2_3;
%...
%... model parameters
%...
%... experiment G3 (catalyst pretreatment first)
     Ptot = 1e5;                    %... total pressure (Pa)
     Tin = 160 + 273.16;            %... inlet temperature (K)
     T0  = 160 + 273.16;            %... initial temperature (K)
     Tw = 29 + 273.16;              %... wall temperature  (K)
     xBin = 0;                      %... mole fraction benzene
     R_kcal	= 1.987;                %... gas constant (kcal/kmole K)
     R_J	= 8.314e3;              %... gas constant (J/kmole K)
     cBin = xBin*Ptot/(R_J*Tin);
     xHin = 1;                      %... mole fraction hydrogen
     cHin = xHin*Ptot/(R_J*Tin);
     xTin = 0;                      %... mole fraction thiophene
     cTin = xTin*Ptot/(R_J*Tin);
     flow = 1551.5e-6/60;           %... total flow rate (m^3/s)
     theta0 = 1;
%...     
     DB = 4e-5;                     %... diffusion coefficients (m^2/s)
     DT = 4e-5;
     eps = 0.6;                     %... void fraction
%...
     MW = 2.106*xHin + 78.12*xBin;                     %... molar weight of the gas mixture (kg/kmole)              
     rhog = MW*273.16*Ptot*0.0075/(22.41*760*Tin);     %... gas density (kg/m^3)  
     cpg = 6.935*xHin + 23.15*xBin;                    %... heat capacity of the gas (kcal/kg K)
     rhocp = 175;                                      %... average heat capacity (kcal/ m^3 K)
     lH2 = 733*0.019/(100*3600);                       %... thermal conductivity of the gas (kcal/s m K)
     rhoc = 354*2/3;                                   %... catalyst density (kg/m^3) with a dilution factor of 2
     Dp = 1.25e-3;                                     %... particle diameter (m) (average between 0.75 and 1.8 mm)
     Rr = 0.822e-2;                                    %... reactor radius (m)
     SA = pi*Rr^2;                                     %... cross-section area (m^2)
     v = flow/(eps*SA);                                %... intersticial velocity (m/s)
     leff = 7*lH2 + 0.8*rhog*cpg*v*eps*Dp;             %... effective thermal conductivity (kcal/s m K)
     alpha = 2.6e-3;                                   %... heat transfer coefficient (kcal/m^2 K)
%...     
%... benzene hydrogenation kinetics
     EB = 13770;                    %... activation energy (kcal/kmole)
     ki0 = 4.22*0.0075;             %... rate constant (kmole/kg s Pa)
     Q = -16470;                    %... heat of adsorption (kcal/kmole)
     K0	= 4.22e-11*0.0075;          %... adsorption constant (1/torr)
%...
%... thiophene poisoning kinetics
     Ed = 1080;                     %... activation energy (kcal/kmole)
     kd0 = (2.40e-2)*0.0075;        %... pre-exponential factor (1/Pa s)
     MT = 1.03e-3;                  %... catalyst adsorption capacity for thiophene (kmole/kg)
%...     
     DH = 49035;                    %... heat of reaction (kcal/kmole)
%...
%... spatial grid
     L = 0.5;                       %... reactor length (m)
     z01 = 0.0;
     zL1 = 0.14;                    %... entrance (inert) section
     z02 = zL1;
     zL2 = z02 + 0.095;             %... catalyst section
     z03 = zL2;
     zL3 = L;                       %... exit (inert) section
     n1 = 71;
     n2 = 71;
     n3 = 71;
     dz1 = (zL1-z01)/(n1-1);
     dz2 = (zL2-z02)/(n2-1);
     dz3 = (zL3-z03)/(n3-1);
     z1 = [z01:dz1:zL1]';
     z2 = [z02:dz2:zL2]';
     z3 = [z03:dz3:zL3]';
     z = [z1;z2;z3];
%...
%... differentiation matrix in zone 1 (convective and diffusive terms)
      D1_1 = five_point_biased_upwind_D1(z1,v);
      D2_1 = five_point_centered_D2(z1);
%...
%... differentiation matrix in zone 2 (convective and diffusive terms)
      D1_2 = five_point_biased_upwind_D1(z2,v);
      D2_2 = five_point_centered_D2(z2);
%...
%... differentiation matrix in zone 3 (convective and diffusive terms)
      D1_3 = five_point_biased_upwind_D1(z3,v);
      D2_3 = five_point_centered_D2(z3);
%...
%... initial conditions
     cB1 = cBin*ones(size(z1));
     cT1 = cTin*ones(size(z1));
     T1 = T0*ones(size(z1));
     cB2 = cBin*ones(size(z2));
     cT2 = cTin*ones(size(z2));
     T2 = T0*ones(size(z2));
     theta = theta0*ones(size(z2));
     cB3 = cBin*ones(size(z3));
     cT3 = cTin*ones(size(z3));
     T3 = T0*ones(size(z3));
%...
%... call to ODE solver
%
%... initial conditions
     x(1:3:3*n1-2,1) = cB1; 
     x(2:3:3*n1-1,1) = cT1; 
     x(3:3:3*n1,1) = T1;
     x(3*n1+1:4:3*n1+4*n2-3,1) = cB2; 
     x(3*n1+2:4:3*n1+4*n2-2,1) = cT2; 
     x(3*n1+3:4:3*n1+4*n2-1,1) = T2;
     x(3*n1+4:4:3*n1+4*n2,1) = theta;
     x(3*n1+4*n2+1:3:3*n1+4*n2+3*n3-2,1) = cB3; 
     x(3*n1+4*n2+2:3:3*n1+4*n2+3*n3-1,1) = cT3; 
     x(3*n1+4*n2+3:3:3*n1+4*n2+3*n3,1) = T3;
%...
     t = 60 * [0:1:30]; % sec
%...
     M = mass;
     options = odeset('Mass',M,'MassSingular','yes','RelTol',1e-3,'AbsTol',1e-3);
%...     
     [tout, yout] = ode15s(@three_zone_reactor_pdes,t,x,options);
%...
%... plot results
     figure(1)
     plot(z1/L,100*yout(:,1:3:3*n1-2)./(yout(:,1:3:3*n1-2)+yout(:,2:3:3*n1-1)+cHin),'k');
     hold on
     plot(z2/L,100*yout(:,3*n1+1:4:3*n1+4*n2-3)./(yout(:,3*n1+1:4:3*n1+4*n2-3)+yout(:,3*n1+2:4:3*n1+4*n2-2)+cHin),'k');
     plot(z3/L,100*yout(:,3*n1+4*n2+1:3:3*n1+4*n2+3*n3-2)./(yout(:,3*n1+4*n2+1:3:3*n1+4*n2+3*n3-2)+yout(:,3*n1+4*n2+2:3:3*n1+4*n2+3*n3-1)+cHin),'k');
     xlabel('z');
     ylabel('xB(z,t)');
     title('Benzene concentration profiles inside the reactor')
%...
     figure(2)
     plot(z1/L,1e4*yout(:,2:3:3*n1-1)./(yout(:,1:3:3*n1-2)+yout(:,2:3:3*n1-1)+cHin),'k');
     hold on
     plot(z2/L,1e4*yout(:,3*n1+2:4:3*n1+4*n2-2)./(yout(:,3*n1+1:4:3*n1+4*n2-3)+yout(:,3*n1+2:4:3*n1+4*n2-2)+cHin),'k');
     plot(z3/L,1e4*yout(:,3*n1+4*n2+2:3:3*n1+4*n2+3*n3-1)./(yout(:,3*n1+4*n2+1:3:3*n1+4*n2+3*n3-2)+yout(:,3*n1+4*n2+2:3:3*n1+4*n2+3*n3-1)+cHin),'k');
     xlabel('z');
     ylabel('xT(z,t)');
     title('Thiophene concentration profiles inside the reactor')
%...
     figure(3)
%     record = avifile('temperature1.avi')
     for k = 1:length(tout),
%       hold off
       plot(z1/L,yout(k,3:3:3*n1)-273.16,'r');
       hold on
       plot(z2/L,yout(k,3*n1+3:4:3*n1+4*n2-1)-273.16,'r');
       plot(z3/L,yout(k,3*n1+4*n2+3:3:3*n1+4*n2+3*n3)-273.16,'r');
       xlabel('z');
       ylabel('T(z,t)');
       title('Temperature profiles inside the reactor')
       axis([0 1 20 180]);
%       mov(k) = getframe(3)
%       record = addframe(record,mov(k))
     end
%     record = close(record)
%...
     figure(4)
     plot(z2/L,yout(:,3*n1+4:4:3*n1+4*n2),'k');
     xlabel('z');
     ylabel('theta(z,t)');
     title('Catalyst activity')
%...
%... read the stopwatch timer
     tcpu=toc;

  
