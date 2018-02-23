%... The MatMol Group (2016)
%... Shock-tube problem  + KT slope limiter
%...
%... u  = -v + eps*u
%...  t     z        zz
%...
%... v  = -[(gam-1)*w-.5*(gam-3)*(v^2)/u]  + eps*v
%...  t                                  z        zz
%...
%... w  = -{[gam*w-.5*(gam-1)*(v^2)/u]*v/u}  + eps*w
%...  t                                    z        zz
%...
%... with   0 < z < 1       t > 0
%...
%...        gam = 1.4 (perfect gas)
%...
%...
%... I.C.:  u(z,0) = 1          0 <= z <= .5-5*eps
%...                   
%...                 linear     .5-5*eps <= z <= .5+5*eps
%...
%...                 0.125      .5+5*eps <= z <= 1
%...
%...        v(z,0) = 0
%...
%...        w(z,0) = 2.5        0 <= z <= .5-5*eps
%...               
%...                 linear     .5-5*eps <= z <= .5+5*eps
%...
%...                 0.25       .5+5*eps <= z <= 1
%...
%...
%... B.C.:  u (0,t) = 0      u (1,t) = 0
%...         z                z
%...
%...        v(0,t) = 0       v(1,t) = 0
%...
%...        w (0,t) = 0      w (1,t) = 0
%...         z                z
%...
     close all
     clear all
%...
%... start a stopwatch timer
     tic
     global n dz z D2 eps gam npdes
%...
%... spatial grid
%...
     z0 = 0;
     zL = 1;
     n = 801;
     nzout = n;
     dz = (zL - z0)/(n-1);
     z = (z0:dz:zL)';
%...
%... problem constants
%...
     eps = .001;
     gam = 1.4;
     npdes = 3;
%...
%... initial conditions
%...
     u = zeros(n,1);
     v = zeros(n,1);
     w = zeros(n,1);
     for i = 1:n
         if z(i) <= .5 - 5*eps
             u(i,1) = 1;
             w(i,1) = 2.5;
         elseif z(i) > .5 - 5*eps && z(i) <= .5 + 5*eps
             u(i,1) = 1 + (.125-1)*(z(i) - .5 + 5*eps)/(10*eps);
             w(i,1) = 2.5 + (.25-2.5)*(z(i) - .5 + 5*eps)/(10*eps);
         else
             u(i,1) = .125;
             w(i,1) = .25;
         end
     end
     x = [u;v;w];
%...
%... differentiation matrix (diffusive term - direct differentiation)
%...
     D2 = five_point_centered_D2(z);
%...  
%... call to ODE solver 
%...     
     time = [0 .15 .23 .28 .32];
%...
     options = odeset('RelTol',1e-6,'AbsTol',1e-6);
     options = odeset(options,'JPattern',Jpattern(n,npdes));
     [tout,yout] = ode15s(@shock_tube_pdes,time,x,options);
%...       
     u = yout(:,1:n);
     v = yout(:,n+1:2*n);
     w = yout(:,2*n+1:3*n);
%...
%... B.C.
%...
     u(:,1) = u(:,2);
     u(:,n) = u(:,n-1);
     v(:,1) = 0;
     v(:,n) = 0;
     w(:,1) = w(:,2);
     w(:,n) = w(:,n-1);
            
     figure(1)
     plot(z,u);
     xlabel('z')
     ylabel('\rho(z,t)')
     axis([0 1 0 1.1])
%...     
     figure(2)
     plot(z,v);
     xlabel('z')
     ylabel('m(z,t)')
%...     
     figure(3)
     plot(z,w);
     xlabel('z')
     ylabel('e(z,t)')
%...     
%... read the stopwatch timer
%... 
     tcpu = toc
