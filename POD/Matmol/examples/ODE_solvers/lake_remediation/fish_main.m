%... The MatMol Group (2016)
%...
%... Effect of liming a fish population in an acidified lake
%...
%... consider a mathematical model to study the effect of acid rains
%... on the fish population of a lake and the effect of remedial liming,
%... which is described in (Ghosh, 2002)
%...
%... the fish population N(t) is growing logistically
%...
%... N  = r*N - r0*N^2/K - H
%...  t  
%...
%... where r(C) is the specific growth rate, which depends on the
%... acid concentration in the following way
%...
%... r(C) = r0                 if C < Clim
%... r(C) = r0-alpha*(C-Clim)  if Clim < C < Cdeath
%... r(C) = 0                  if Cdeath < C < Q/delta
%...
%... K(C) is the carrying capacity, which also depends on C
%...
%... K(C) = K0                 if C < Clim
%... K(C) = K0-beta*(C-Clim)   if Clim < C < Cdeath
%... K(C) = Klim               if Cdeath < C < Q/delta
%... 
%... and H is the harvesting rate
%...
%... Clim denotes the critical value of the acid concentration (between 0 and 
%... Clim, acid is harmless to the fish population)
%...
%... the acid concentration is described by
%...
%... C  = Q - delta*C - delta0*E
%...  t  
%... 
%... where Q is the inlet acid flow rate (due to acid rains), delta is the 
%... natural depletion rate, whereas delta0*E represents the effect of liming
%...
%... E is the liming effort applied to maintain the lake at a permissible
%... acid concentration Clim
%...
%... E  = eta*(C-Clim) - eta0*E
%...  t  
%...
%... where eta0 is the natural depletion rate of E
%...
%... the initial conditions (ICs) are taken as
%...
%... N(t=0) = 72500    C(t=0) = 80    E(t=0) = 190          
%...
%... Reference
%...
%... M. Ghosh, 
%... Effect of liming a fish population in an acidified lake:
%... a simple mathematical model, 
%... Applied Mathematics and Computation (2002), in press.
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
     global r0 Clim alpha Cdeath K0 Klim beta H Q delta delta0 eta eta0
     global fac thresh vectorized
%...
%... model parameters
     r0 = 0.02;
     Clim = 50;
     alpha = 0.0001;
     Cdeath = (r0+alpha*Clim)/alpha;
     K0 = 100000;
     Klim = 100;
     beta = 0.05; 
     H = 100;
     Q = 2;
     delta = 0.002;
     delta0 = 0.005;
     eta0 = 0.004;
%...  
%... select the control action 'eta' (comment/decomment one of the actions)
%...
%     action = 'strong'         
%...
%      action = 'moderate'           
%...
      action = 'weak'     
%...
     switch action    
%...
%... strong
       case('strong')
           eta = 0.5;
%... moderate
       case('moderate')
           eta = 0.1;
%... weak
       case('weak')
           eta = 0.04;     
%...
     end
%...
%... equilibrium points
     [N1star,N2star] = equilibrium(r0,Clim,alpha,Cdeath,K0,Klim,beta,H,Q,delta,delta0,eta,eta0)
%...
%... initial conditions
     t0=0;
     tf=3560;
     N = 72500;
%     N = 5000;
     C = 80;
     E = 190;
     x=[N C E]';
%...  
%... call to ODE solver (comment/decomment one of the methods to select a solver)
%...
%      method = 'euler'         
%...
%      method = 'midpoint'         
%...
%      method = 'heun'           
%...
%      method = 'heun12'           
%...
%      method = 'rk4'     
%...
%      method = 'rkf45'     
%...
%      method = 'ros3p'     
%...
%      method = 'ros23p'     
%...
      method = 'ode45'     
%...
%      method = 'ode15s'
%...
%      method = 'lsodes'
%...
     switch method    
%...
%... Euler
       case('euler')
           Dt=0.1;
           Dtplot=20;
           [tout, xout] = euler_solver(@fish_odes, t0, tf, x, Dt, Dtplot);
%... midpoint
       case('midpoint')
           Dt=5;
           Dtplot=20;
           [tout, xout] = midpoint_solver(@fish_odes, t0, tf, x, Dt, Dtplot);
%... Heun
       case('heun')
           Dt=5;
           Dtplot=20;
           [tout, xout] = heun_solver(@fish_odes, t0, tf, x, Dt, Dtplot);
%... Heun12
       case('heun12')
           hmin=0.0001;
           nstepsmax=1e5;
           abstol=1e-3;
           reltol=1e-3;
           Dtplot=20;
           [tout, xout] = heun12_solver(@fish_odes,t0,tf,x,hmin,nstepsmax,abstol,reltol,Dtplot);
%... rk4
       case('rk4')
           Dt=5;
           Dtplot=20;
           [tout, xout] = rk4_solver(@fish_odes, t0, tf, x, Dt, Dtplot);
%... rkf45
       case('rkf45')
           hmin=0.0001;
           nstepsmax=1e5;
           abstol=1e-3;
           reltol=1e-3;
           Dtplot=20;
           [tout, xout, eout] = rkf45_solver(@fish_odes,t0,tf,x,hmin,nstepsmax,abstol,reltol,Dtplot);
           figure(2)
           plot(tout/356,eout(:,1),':r');
           hold on
           plot(tout/356,eout(:,2),'--b');
           plot(tout/356,eout(:,3),'--g');
           xlabel('t');
           ylabel('e(t)');
%... ros3p
       case('ros3p')
           Dt=1;
           Dtplot=20;
           fac = [];
           thresh = 1e-12;
%           [tout, xout] = ros3p_solver(@fish_odes,@jacobian,t0,tf,x,Dt,Dtplot);
           [tout, xout] = ros3p_solver(@fish_odes,@jacobian_num,t0,tf,x,Dt,Dtplot);
%... ros23p
       case('ros23p')
           hmin=0.0001;
           nstepsmax=1e5;
           abstol=1e-3;
           reltol=1e-3;
           Dtplot=1;
           fac = [];
           thresh = 1e-12;
%           [tout, xout, eout] = ros23p_solver(@fish_odes,@jacobian,@ft,t0,tf,x,hmin,nstepsmax,abstol,reltol,Dtplot);
           [tout, xout, eout] = ros23p_solver(@fish_odes,@jacobian_num,@ft,t0,tf,x,hmin,nstepsmax,abstol,reltol,Dtplot);
           figure(2)
           plot(tout,eout(:,1),':r');
           hold on
           plot(tout,eout(:,2),'--b');
           plot(tout,eout(:,3),'--g');
           xlabel('t');
           ylabel('e(t)');
%... ode45
       case('ode45')
           options = odeset('RelTol',1e-3,'AbsTol',1e-3);
           Dtplot=20;
           t=[t0:Dtplot:tf];
           [tout, xout] = ode45(@fish_odes,t,x,options);
%... ode15s
       case('ode15s')
           options = odeset('RelTol',1e-3,'AbsTol',1e-3);
           Dtplot=20;
           t=[t0:Dtplot:tf];
           [tout, xout] = ode15s(@fish_odes,t,x,options);
%... lsodes
       case('lsodes')
            resname  = 'fish_odes'; 
            jacname  = '[]';
            neq = 3;
            Dtplot=20;
            tlist = [t0+Dtplot:Dtplot:tf];
            itol   = 1; 
            abstol = 1e-3;
            reltol = 1e-3;
            itask  = 1; 
            istate = 1;
            iopt   = 0; 
            lrw = 50000; 
            rwork = zeros(lrw,1); 
            liw = 50000; 
            iwork = zeros(liw,1); 
            mf = 222; 
            [tout, xout] = Lsodes(resname,jacname,neq,x,t0,tlist,itol,reltol,abstol,itask,istate,iopt,rwork,lrw,iwork,liw,mf); 
%...
     end
%...
%... plot results
     figure(1)
     subplot(3,1,1)
     plot(tout/356,xout(:,1),'-');
%     xlabel('t [years]');
     ylabel('N(t)');
%     title('Fish population');
     subplot(3,1,2)
     plot(tout/356,xout(:,2),'-');
%     xlabel('t [years]');
     ylabel('C(t)');
%     title('Acid concentration')
     subplot(3,1,3)
     plot(tout/356,xout(:,3),'-');
     xlabel('t [years]');
     ylabel('E(t)');
%     title('Liming effort');
%...
%... read the stopwatch timer
     tcpu=toc;

  
