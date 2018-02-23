%... The MatMol Group (2016)
%...
%... Brusselator
%...
%...           2
%... u  = A + u v - (B+1)u + alpha*u
%...  t                             zz
%...
%...            2
%... v  = Bu - u v + alpha*v
%...  t                     zz
%...
%...
%... where
%...
%...    A = 1
%...
%...    B = 3
%...  
%...    alpha = 1/50
%...
%...    0 < z <1
%...
%...    t > 0
%...
%... 
%... with
%...
%...    C.L. :  u(0,t) = u(1,t) = 1
%...
%...            v(0,t) = v(1,t) = 3
%...
%...    C.I. :  u(z,0) = 1 + sin(2*pi*z)
%...
%...            v(z,0) = 3
     clear all
     close all
%...
%... set global variables
     global A B alpha n dz ul ur vl vr 
%...
%... model parameters
     A = 1;
     B = 3;
     alpha = 1/50;
     ul = 1;
     ur = 1;
     vl = 3;
     vr = 3;
%...
%... spatial grid
     n = 101;
     dz = 1/(n-1);
     z = [0:dz:1]';
%...
%... initial conditions
     u = 1 + sin(2*pi*z);
     v = 3*ones(n,1);
     x0 = [u ; v];
%...  
%... Call to ODE solver
     t0 = 0;
     tf = 20;
     fac = [];
     thresh = 1e-16;
%...
    Dtplot = 0.5;
    options = odeset('RelTol',1e-3,'AbsTol',1e-3);
    options = odeset(options,'Jacobian',@jacobian_num);
%    options = odeset(options,'Jacobian',@jacobian_complex_diff);
    S = jpattern_num;
%     S = jpattern_complex_diff;
     options = odeset(options,'JPattern',S)
%...
    [tout,xout] = ode15s(@brusselator_pdes,[t0:Dtplot:tf],x0,options);
%...
%... separate u and v from x
     uout = xout(:,1:n);
     vout = xout(:,n+1:2*n);
%...               
%... plot results
     figure
     mesh(z,tout,uout);
     hold;
     axis([0 1 0 20 0 3]);
     xlabel('z');
     ylabel('t');
     zlabel('u(z,t)');
     title('Brusselator')
%...     
     figure
     mesh(z,tout',vout);
     hold;
     axis([0 1 0 20 1 4.5]);
     xlabel('z');
     ylabel('t');
     zlabel('v(z,t)');
     title('Brusselator')



