%... The MatMol Group (2016)
     function xt = three_zone_reactor_pdes(t,x)
%...
%... set global variables
     global v DB DT rhocp cpg rhog Tw Rr DH lH2 Dp leff eps alpha rhoc;
     global ki0 K0 Q EB R_kcal R_J kd0 Ed MT;
     global Ptot cBin cTin Tin cHin xHin;
     global z01 zL1 z02 zL2 z03 zL3 n1 n2 n3 z1 z2 z3 D1_1 D2_1 D1_2 D2_2 D1_3 D2_3;
%...
%... Transfer dependent variables  
     cB1 = x(1:3:3*n1-2,1);
     cT1 = x(2:3:3*n1-1,1);
     T1 = x(3:3:3*n1,1);
%...     
     cB2 = x(3*n1+1:4:3*n1+4*n2-3,1);
     cT2 = x(3*n1+2:4:3*n1+4*n2-2,1);
     T2 = x(3*n1+3:4:3*n1+4*n2-1,1);
     theta = x(3*n1+4:4:3*n1+4*n2,1);
%...     
     cB3 = x(3*n1+4*n2+1:3:3*n1+4*n2+3*n3-2,1);
     cT3 = x(3*n1+4*n2+2:3:3*n1+4*n2+3*n3-1,1);
     T3 = x(3*n1+4*n2+3:3:3*n1+4*n2+3*n3,1);
%...
%... spatial derivatives - 1st zone
     cB1z = D1_1*cB1;
     cT1z = D1_1*cT1;
     T1z = D1_1*T1;
%...
     cB1zz = D2_1*cB1;
     cT1zz = D2_1*cT1;
     T1zz = D2_1*T1;
%...
%... spatial derivatives - 2nd zone
     cB2z = D1_2*cB2;
     cT2z = D1_2*cT2;
     T2z = D1_2*T2;
%...
     cB2zz = D2_2*cB2;
     cT2zz = D2_2*cT2;
     T2zz = D2_2*T2;
%...
%... spatial derivatives - 3rd zone
     cB3z = D1_3*cB3;
     cT3z = D1_3*cT3;
     T3z = D1_3*T3;
%...
     cB3zz = D2_3*cB3;
     cT3zz = D2_3*cT3;
     T3zz = D2_3*T3;
%...
%... several operating conditions
%...
%... 1) catalyst pretreatment with hydrogen at 160°C (for 20 min)
%...
     rB = 0;                           %... reaction rates in 2nd zone
     rd = 0;
     rT = 0;
%...
%... 2) benzene hydrogenation (experiment G3)
%...
     if (t > 20*60) & (t < 30*60)
%...         
     Tin = 160 + 273.16;              %... inlet temperature (K)
     xBin = 2*0.033;                   %... mole fraction benzene
     cBin = xBin*Ptot/(R_J*Tin);
     xHin = 1-xBin;                    %... mole fraction hydrogen
     cHin = xHin*Ptot/(R_J*Tin);
     xTin = 0;                         %... mole fraction thiophene
     cTin = xTin*Ptot/(R_J*Tin);
%...
     MW = 2.106*xHin + 78.12*xBin;                     %... molar weight of the gas mixture (kg/kmole)              
     rhog = MW*273.16*Ptot*0.0075/(22.41*760*Tin);     %... gas density (kg/m^3)  
     cpg = 6.935*xHin + 23.15*xBin;                    %... heat capacity of the gas (kcal/kg K)
     leff = 7*lH2 + 0.8*rhog*cpg*v*eps*Dp;             %... effective thermal conductivity (kcal/s m K)
%...
     xB2 = cB2./(cB2+cT2+cHin);                        %... reaction rates in 2nd zone
     rB = ki0*K0*Ptot^2*(xHin*xB2.*theta.*exp((-Q-EB)./(R_kcal*T2)))./(1+K0*Ptot*xB2.*exp(-Q./(R_kcal*T2)));
     rd = 0;
     rT = 0;
% %...
% %... 3) catalyst poisoning (experiment G3)
% %...
%      elseif t > 60*120
% %...         
%      Tin = 160 + 273.16;              %... inlet temperature (K)
%      xBin = 2*0.033;                   %... mole fraction benzene
%      cBin = xBin*Ptot/(R_J*Tin);
%      xHin = 1-xBin;                    %... mole fraction hydrogen
%      cHin = xHin*Ptot/(R_J*Tin);
%      xTin = 1.136*xBin/100;            %... mole fraction thiophene
%      cTin = xTin*Ptot/(R_J*Tin);
% %...
%      MW = 2.106*xHin + 78.12*xBin;                     %... molar weight of the gas mixture (kg/kmole)              
%      rhog = MW*273.16*Ptot*0.0075/(22.41*760*Tin);     %... gas density (kg/m^3)  
%      cpg = 6.935*xHin + 23.15*xBin;                    %... heat capacity of the gas (kcal/kg K)
%      leff = 7*lH2 + 0.8*rhog*cpg*v*eps*Dp              %... effective thermal conductivity (kcal/s m K)
% %...
%      xB2 = cB2./(cB2+cT2+cHin);                        %... reaction rates in 2nd zone
%      xT2 = cT2./(cB2+cT2+cHin);
%      rB = ki0*K0*Ptot^2*(xHin*xB2.*theta.*exp((-Q-EB)./(R_kcal*T2)))./(1+K0*Ptot*xB2.*exp(-Q./(R_kcal*T2)));
%      rd = kd0*Ptot*xT2.*theta.*exp(-Ed./(R_kcal*T2));
%      rT = MT*rd;
% %...     
     end
%...
%... temporal derivatives
%...
     cB1t = -v*cB1z + DB*cB1zz;
     cT1t = -v*cT1z + DT*cT1zz;
     T1t = -((eps*v*rhog*cpg)/rhocp)*T1z + (leff/rhocp)*T1zz + 2*alpha*(Tw - T1)/(Rr*rhocp);
%...     
     cB2t = -v*cB2z + DB*cB2zz - rhoc*rB/eps;
     cT2t = -v*cT2z + DT*cT2zz - rhoc*rT/eps;
     T2t = -((eps*v*rhog*cpg)/rhocp)*T2z + (leff/rhocp)*T2zz + 2*alpha*(Tw - T2)/(Rr*rhocp) + (DH/rhocp)*rhoc*rB;
     thetat = -rd;
%...     
     cB3t = -v*cB3z + DB*cB3zz;
     cT3t = -v*cT3z + DT*cT3zz;
     T3t = -((eps*v*rhog*cpg)/rhocp)*T3z + (leff/rhocp)*T3zz + 2*alpha*(Tw - T3)/(Rr*rhocp);
%...
%... boundary conditions at z = z01
     cB1t(1) = cBin - cB1(1);
     cT1t(1) = cTin - cT1(1);
     T1t(1) = Tin - T1(1);
%...
%... boundary conditions at z = zL1 = z02
     cB1t(n1) = cB2z(1) - cB1z(n1);
     cT1t(n1) = cT2z(1) - cT1z(n1);
     T1t(n1) = T2z(1) - T1z(n1);
%...
     cB2t(1) = cB1(n1) - cB2(1);
     cT2t(1) = cT1(n1) - cT2(1);
     T2t(1) = T1(n1) - T2(1);
%...
%... boundary conditions at z = zL2 = z03
     cB2t(n2) = cB3z(1) - cB2z(n2);
     cT2t(n2) = cT3z(1) - cT2z(n2);
     T2t(n2) = T3z(1) - T2z(n2);
%...
     cB3t(1) = cB2(n2) - cB3(1);
     cT3t(1) = cT2(n2) - cT3(1);
     T3t(1) = T2(n2) - T3(1);
%...
%... boundary conditions at z = zL
     cB3t(n3) = - cB3z(n3);
     cT3t(n3) = - cT3z(n3);
     T3t(n3) = - T3z(n3);
%...
%... Transfer temporal derivatives
     xt(1:3:3*n1-2,1) = cB1t; 
     xt(2:3:3*n1-1,1) = cT1t; 
     xt(3:3:3*n1,1) = T1t;
%...     
     xt(3*n1+1:4:3*n1+4*n2-3,1) = cB2t; 
     xt(3*n1+2:4:3*n1+4*n2-2,1) = cT2t; 
     xt(3*n1+3:4:3*n1+4*n2-1,1) = T2t;
     xt(3*n1+4:4:3*n1+4*n2,1) = thetat;
%...     
     xt(3*n1+4*n2+1:3:3*n1+4*n2+3*n3-2,1) = cB3t; 
     xt(3*n1+4*n2+2:3:3*n1+4*n2+3*n3-1,1) = cT3t; 
     xt(3*n1+4*n2+3:3:3*n1+4*n2+3*n3,1) = T3t;
