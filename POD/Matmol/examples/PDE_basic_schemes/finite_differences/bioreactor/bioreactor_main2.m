%... The MatMol Group (2016)
%... Dynamic analysis of a bioreactor (block storage)
%...
%... The following equations model a bioreactor
%...
%...   Mass balances
%...
%...   S  = -v*S   - nu1*(1-eps)*phi1/eps                          (1)
%...    t       z
%...
%...
%...   X  = phi1 - phi2                                            (2)
%...    t
%...
%...
%... The variables and parameters for this model are
%...
%...   Substrate                       S
%...
%...   Biomass                         X
%...
%...   Time                            t
%...
%...   Axial position                  z
%...
%...   Bioreactor length               zL
%...
%...
     close all
     clear all
%...
%... start a stopwatch timer
     tic
%...
%... set global variables
     global v nu1 eps mumax Kc kd;
     global Sin gamma
     global z0 zL z nz D1z;
%...
%... model parameters
     eps = 0.5;
     A = 0.04;
     F = 2e-3;
     v = F/(eps*A);
     nu1 = 0.4;
     mumax = 0.35;
     Kc = 0.4;
     kd = 0.05;
%...
%... inlet concentration
     Sin = 5;
     gamma = 10;
%...
%... initial conditions
     S0 = 0;
     X0 = 10;
%...
%... grid in axial direction
     z0 = 0.0;
     zL = 1.0;
     nz = 101;
     dz = (zL-z0)/(nz-1);
     z = [z0:dz:zL]';
%...
%... initial conditions
     S = S0*ones(nz,1);
     X = X0*ones(nz,1);
%...
%... transfer dependent variables
     x(1:nz) = S; 
     x(nz+1:2*nz) = X;
%...
%... differentiation matrix in z (convective term)
     D1z = two_point_upwind_uni_D1(z0,zL,nz,A);
%...  
%... call to ODE solver 
%...
     t0 = 0.0;
     tf = 50.0;
     Dtplot = 1.0;
     options = [];
%...     
     options = odeset('RelTol',1e-3,'AbsTol',1e-3);
%...
%... compute the Jacobian matrix
%     options = odeset(options,'Jacobian',@jacobian_num2);
%     options = odeset(options,'Jacobian',@jacobian_complex_diff2);
%...
%... or compute the sparcity pattern
%     Sp = jpattern_num2;
%     Sp = jpattern_complex_diff2;
%     options = odeset(options,'JPattern',Sp)
%...     
     [tout, yout] = ode15s(@bioreactor_pdes2,[t0:Dtplot:tf],x,options);
%...     
     for k=1:length(tout),
         S_out=yout(k,1:nz);
         figure(1);
         plot(z,S_out);
         xlabel('z');
         ylabel('S(z,t)');
         title('Substrate')
         axis([0 1 0 7]);
%         pause(.1)
     end
    display('Biomass display starts soon'); pause(1)
     for k=1:length(tout),
         X_out=yout(k,nz+1:2*nz);
         figure(2);
         plot(z,X_out);
         xlabel('z');
         ylabel('X(z,t)');
         title('Biomass')
         axis([0 1 0 70]);
%         pause(.1)
     end     
%...
%... read the stopwatch timer
     tcpu=toc;
