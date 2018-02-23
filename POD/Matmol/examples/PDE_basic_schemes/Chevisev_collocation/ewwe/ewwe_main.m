%...  The Matmol group (2016)
%...
%... Generalized equal-width wave equation (implicit solution with ode15s)
%...
%...
%... x  + a*(x^p)*x  - mu*x    = 0
%...  t            z       zzt
%...
%...
%... where
%...
%...  x       dependent variable
%...
%...  t       time
%...
%...  z       space
%...
%...  mu      positive real constant
%...
%...  p       1 for the classical EWWE
%...
%... Reference
%...
%... S. Hamdi, W.H. Enright, W.E. Schiesser, J.J. Gottlieb
%... Exact solutions and invariants of motion for general types of
%... regularized long wave equations
%... Mathematics and Computers in Simulation 65 (2004), 535-545
%...
     close all
     clear all
%...
%... start a stopwatch timer
     tic
%...
%... set global variables
     global a p mu c zi;
     global z0 zL z n D1;
%...
%... spatial grid
     z0 = -100.0;
     zL = 100.0;
     n = 200;       % spectral method
%     n = 500;       % FDs
%...
%... spectral differentiation
     [D1,z] = periodic_spectral_D1(z0,zL,n);
%... possibly compared to FDs
%     D1=three_point_centered_uni_D1(z0,zL,n);
%...
%... parameters and initial conditions
     p = 1;
     a = 1.5;
     mu = 3;
     c = 0.6;
     zi = 0;
     c2 = 0.1;
     z01 = 0;
     z02 = 30;
%...
%... one solitary wave
     x = (((p+1)*(p+2)*c*(sech(0.5*p*sqrt(1/mu)*z)).^2)/(2*a)).^(1/p);
%...
%... interaction between two solitary waves
%     x = (((p+1)*(p+2)*c*(sech(0.5*p*sqrt(1/mu)*(z-z01))).^2)/(2*a)).^(1/p) +...
%         (((p+1)*(p+2)*c2*(sech(0.5*p*sqrt(1/mu)*(z-z02))).^2)/(2*a)).^(1/p);    
%...  
%... call to ODE solver 
%...
%... ode15s
     t0 = 0;
     dt = 10;
     tf = 100;
     t = [t0:dt:tf];
     options = odeset('RelTol',1e-6,'AbsTol',1e-8);
     M = mass;
     options = odeset(options,'Mass',M,'MassSingular','no');
%     options = odeset(options,'JPattern',jpattern(n));
%...
%     figure(1);
%     colordef black;
%     ind = [1 round(n/4) round(n/2) round(3*n/4) n];
%     options = odeset(options,'OutputFcn',@odeplot,'OutputSel',[ind]);
%...     
     [tout, xout] = ode15s(@ewwe_pde,t,x,options);
%...
%... computation of solution invariants
     for k=1:length(tout),
%...         
         I1(k) = trapz(z,xout(k,:));
         xoutz = D1*xout(k,:)';
         I2(k) = trapz(z,xout(k,:)'.^2 + mu*xoutz.*xoutz);
         I3(k) = trapz(z,(2*a/(p+2))*xout(k,:).^(p+2));
%...
%...      exact values (for z0 = 0, p = 1, a = 1.5 and for a finite, symetric [-L +L], spatial domain)
%...      
%          I1_exact(k) = 12*c*sinh(0.5*sqrt(1/mu)*(zL-z0)/2)/( a*sqrt(1/mu)*cosh(0.5*sqrt(1/mu)*(zL-z0)/2) );
%          I2_exact(k) = 36/5*c^2*sinh(0.5*sqrt(1/mu)*(zL-z0)/2) * (  2*(cosh(0.5*sqrt(1/mu)*(zL-z0)/2) )^2 ...
%                                                                   + 4*(cosh(0.5*sqrt(1/mu)*(zL-z0)/2) )^4 ...
%                                                                   - 1                                    ) ...
%                         / ( a^2*sqrt(1/mu)*( cosh(0.5*sqrt(1/mu)*(zL-z0)/2))^5 );
%          I3_exact(k) = 36/5*c^3*sinh(0.5*sqrt(1/mu)*(zL-z0)/2) * (  4*(cosh(0.5*sqrt(1/mu)*(zL-z0)/2) )^2 ...
%                                                                   + 8*(cosh(0.5*sqrt(1/mu)*(zL-z0)/2) )^4 ...
%                                                                   + 3                                    ) ...
%                         / ( a^3*sqrt(1/mu)*( cosh(0.5*sqrt(1/mu)*(zL-z0)/2))^5 );
%...
%...      asymptotic values (for z0= 0, for p = 1, a = 1.5 and for an infinite dimensional domain)
          I1_exact(k) = 12*c*sqrt(mu)/a;
          I2_exact(k) = 144/5*c^2*sqrt(mu)/a^2;
          I3_exact(k) = 288/5*c^3*sqrt(mu)/a^3;
     end    
%...
%... plot results
     set(0,'defaultaxesfontsize',12,'defaultaxeslinewidth',.7,'defaultlinelinewidth',1.5,'defaultpatchlinewidth',.7);
%...
     figure(2);
     plot(z,xout,'k');
     xlabel('z');
     ylabel('x(t,z)');
     hold on
%...     
     for k=1:length(tout),
         for i=1:n,
             xexact(i) = ewwe_exact(z(i),tout(k));
         end;
         plot(z,xexact,':r')
     end;
%...
     figure(3)
%     mesh(z,tout,xout);
     surf(z,tout,xout);
     view(10,80);
     grid off;
     axis([z0 zL t0 tf min(min(xout)) max(max(xout))]);
     xlabel('z');
     ylabel('t');
     zlabel('u(t,z)');
%...
     figure(4)
     plot(tout,I1,'k')
     hold on
     plot(tout,I2,'b')
     plot(tout,I3,'g')
     
     plot(tout,I1_exact,'dr')
     plot(tout,I2_exact,'or')
     plot(tout,I3_exact,'+r')
     xlabel('t');
     ylabel('I_1, I_2, I_3');
     axis([0 tf 0 20]);
     legend('I_1','I_2','I_3');
     title('Invariants of Motion')
%...
%... read the stopwatch timer
     tcpu=toc;

  
