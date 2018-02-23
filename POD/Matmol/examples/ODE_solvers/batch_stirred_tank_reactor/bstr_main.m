%... The MatMol Group (2016)
%...
%... Batch Stirred Tank Reactor
%...
%... consider a BSTR in which the following reactions take place
%... (modified from Economou, Morari and Palsson, 1986)
%...
%... A <----> B
%...
%... the reaction rates are linear in either the component A or B,
%... and follow the Arrhenius temperature dependency
%...
%... kd = Cd*exp(-Ed/(R*T))
%...
%... kr = Cr*exp(-Er/(R*T))
%...
%... the material balances are given by
%...
%... A  = - kd*A + kr*B
%...  t  
%...
%... B  = kd*A - kr*B
%...  t  
%...
%... the energy balance is given by
%...
%... T  = (kd*A - kr*B)*DH/(rho*cp)
%...  t   
%...
%... where
%...
%...  A       component concentration
%...
%...  B       component concentration
%...
%...  T       temperature
%...
%...  t       time
%...
%...  V       reactor volume
%...
%...  rho     density
%...
%...  cp      heat capacity
%...
%...  F       outlet flow rate
%...
%...  Cd      frequency constant (direct reaction)
%...
%...  Cr      frequency constant (reverse reaction)
%...
%...  Ed      activation energy (direct reaction)
%...
%...  Er      activation energy (reverse reaction)
%...
%...  R       gas constant
%...
%...
%... the initial conditions (ICs) are taken as
%...
%... A(t=0) = A0 = 1    B(t=0) = B0 = 0    T(t=0) = T0 = 430          
%...
%... Reference
%...
%... C.G. Economou, M. Morari and B.O. Palsson, 
%... Internal Model Control. 5. Extension to Nonlinear Systems, 
%... Ind. Eng. Chem. Process Des. Dev. 25 (1986), 403-411.
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
     global rhocp Cd Cr Ed Er DH R
%...
%... model parameters
     rhocp = 1000;
     Cd = 5000;
     Cr = 1000000;
     Ed = 10000; 
     Er = 15000; 
     DH = 5000;
     R = 1.987;
%...
%... initial conditions
     t0 = 0;
     tf = 100;
     A0 = 1;
     B0 = 0;
     T0 = 430;
     x = [A0 B0 T0]';
%...  
%... call to ODE solver (comment/decomment one of the methods to select a solver)
%...
%      method = 'euler'         
%...
%...   method = 'heun'           
%...
%   method = 'leap-frog'           
%...
      method = 'euler/leap-frog'           
%...
%      method = 'ode45'     
%...
%      method = 'ode15s'     
%...
     switch method    
%...
%... Euler
       case('euler')
           Dt=1;
           Dt=0.5;
           Dtplot=2;
           Dtplot=1;
           [tout, yout] = euler_solver(@bstr_odes, t0, tf, x, Dt, Dtplot);
%... modified Euler
       case('heun')
           Dt=2;
           Dtplot=2;
          [tout, yout] = heun_solver(@bstr_odes, t0, tf, x, Dt, Dtplot);
%... leap-frog
       case('leap-frog')
           Dt=0.05;
           Dtplot=0.05;
           xm1=x;
           [tout, yout] = leap_frog_solver(@bstr_odes, t0, tf, x, xm1, Dt,Dtplot);
%... euler/leap-frog
       case('euler/leap-frog')
           Dt=0.05;
           Dtplot=0.05;
           [t1, y1] = euler_solver(@bstr_odes, t0, t0+Dt, x, Dt, Dt);
           x0       = y1(end,:)';
           [t2, y2] = leap_frog_solver(@bstr_odes, t0+Dt, tf, x0, x, Dt,Dtplot);
           tout = [t1;t2];
           yout = [y1;y2];
%... ode45
       case('ode45')
           options = odeset('RelTol',1e-3,'AbsTol',1e-3);
           Dtplot=10;
           t=[t0:Stplot:tf];
           [tout, yout] = ode45(@bstr_odes,t,x,options);
%... ode15s
       case('ode15s')
           options = odeset('RelTol',1e-3,'AbsTol',1e-3);
           Dtplot=0.075;
           t=[t0:Dplot:tf];
           [tout, yout] = ode15s(@bstr_odes,t,x,options);
%...
     end
%...     
%... plot results
%...
     subplot(3,1,1)
     plot(tout,yout(:,1),'k');
%     xlabel('t');
     ylabel('A(t)');
%     title('Component A concentration');
     subplot(3,1,2)
     plot(tout,yout(:,2),'k');
%     xlabel('t');
     ylabel('B(t)');
%     title('Component B concentration')
     subplot(3,1,3)
     plot(tout,yout(:,3),'k');
     xlabel('t');
     ylabel('T(t)');
%     title('Temperature');
%...
%... read the stopwatch timer
     tcpu=toc;

  
