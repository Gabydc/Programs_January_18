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
     n = 100;
     dz = 1/n;
     z = [0:dz:1]';
%...
%... initial conditions
     u = 1 + sin(2*pi*z);
     v = 3*ones(n+1,1);
     x0(1:2*n-2,1) = [u(2:n) ; v(2:n)];
%...  
%... Call to ODE solver
     t0 = 0;
     tf = 20;
     fac = [];
     thresh = 1e-16;
%...
     hmin = 0.001;
     nstepsmax = 1e5;
     abstol = 1e-3;
     reltol = 1e-3;
     Dtplot = 0.5;
     [tout,xout,eout] = ros23p_solver(@brusselator_pdes2,@jacobian,@ft,t0,tf,x0,hmin,nstepsmax,abstol,reltol,Dtplot);
%...
%... separate u and v from x
     for i =1:length(tout)
         uf(i,2:n) = xout(i,1:n-1);
         vf(i,2:n) = xout(i,n:2*n-2);
     end
     uf(1:length(tout),1) = 1;
     uf(1:length(tout),n+1) = 1;
     vf(1:length(tout),1) = 3;
     vf(1:length(tout),n+1) = 3;
%...               
%... plot results
     figure
     mesh(z,tout,uf);
     hold;
     axis([0 1 0 20 0 3]);
     xlabel('z');
     ylabel('t');
     zlabel('u(z,t)');
     title('Brusselator')
%...     
     figure
     mesh(z,tout',vf);
     hold;
     axis([0 1 0 20 1 4.5]);
     xlabel('z');
     ylabel('t');
     zlabel('v(z,t)');
     title('Brusselator')
