%... The MatMol Group (2016)
%...
%... Burgers equation
%...
%... The well known Burgers equation :
%...
%... x  = -x*x  + mu*x   
%...  t       z       zz
%...
%...
%... where
%...
%...  x      dependent variable
%...
%...  t       time
%...
%...  z       space
%...
%...  mu      diffusion coefficient
%...
%...
%... an exact solution is given by
%...
%... x(z,t) = burgers_exact(z,t)
%...
%... the boundary conditions (BCs) for Burgers' equation are taken as
%...
%... x(z=zL,t) = burgers_exact(zL,t)     (Dirichlet BCs)
%...
%... x(z=zR,t) = burgers_exact(zR,t)               
%...
%... where
%...
%...  zL < z < zR is the spatial domain
%...
%... the initial condition (IC) for Burgers' equations are taken as
%...
%... x(z,t=0) = burgers_exact(z,0)
%...
%... Implementation of MOVGRID - Finite differences version
%...
%...
     close all
     clear all
%...
%... start a stopwatch timer
%...
     tic
%... 
%... set global variables
%...
     global zL zR n ne X eq choice
     global alpha kappa mu tau
     global u z
     global visc
%...
%... temporal conditions
     t=[0:0.1:1];
%...
%... moving grid parameters
     alpha=1;
     kappa=2;
     mu=kappa*(kappa+1);
     tau=1e-4;
%...
%... PDE parameters
     ne=1;
     visc=0.001;
%...
%... initial spatial grid z : if n is the number of moving grid points,
%... dim(z) = n+2
     zL=0.0;
     zR=1.0;
     n=200;
     dz=(zR-zL)/(n+1);
     z=[zL:dz:zR]';
%...
%... initial dependent variables u(n+2,ne)
     for j=1:ne,
         for i=1:n+2,
             u(i,j) = burgers_exact(z(i),0);
         end
     end
%...
%... global variables for JPat : function JPat implements the sparsity
%... pattern of the jacobian in X ; informations depending on the
%... choosen spatial derivative matrices are stocked in the ne eq(ne,ne,3) matrices.
%... 
%...    example : hypothetic problem : 
%...
%...        u1  = -f1  + eps1*u1   + k1*u2      where  f1 depends on u3
%...          t      z          zz
%...
%...        u2  = -f2  + eps2*u2                where  f2 depends on u2 and u3
%...          t      z          zz
%...
%...        u3  = -f3  + eps3*u3                where  f3 depends on u1, u2 and u3
%...          t      z          zz
%...
%...    where the 1st (2nd) derivative is computed by a 3 (5) pts centered finite difference
%... 
%...    It results :
%...
%...        eq(:,:,1) = [-2 2 1; 0 0 1;-1 1 1];
%...        eq(:,:,2) = [ 0 0 0;-2 2 1;-1 1 1];
%...        eq(:,:,3) = [-1 1 1;-1 1 1;-2 2 1];

     X=zeros(n*(ne+1),n*(ne+1));
%... 
%... see right_vect file :
%... diffusion term : D2 = 3 pts centered finite difference
%... convection term : D1 = 5 pts biased_upwind finite difference
     eq(:,:,1) = [-3 1 1];
%...
%... monitoring choice : 
%...
%...                                                     1    ne (u(i+1,j)-u(i,j))^2
%...    choose between : 'der1'  : m(i) = sqrt[ alpha + ---  sum ------------------- ]
%...                                                     ne  j=1   (z(i+1)-z(i))^2
%...
%...        
%...
%...                                                     1    ne     uzz(i+1,j)+uzz(i,j)
%...                     'der2'  : m(i) = sqrt[ alpha + ---  sum abs(-------------------)]
%...                                                     ne  j=1           2  
%...
%...
%...
%...                                                     1    ne (fl[u(i+1,j)]-fl[u(i,j)])^2
%...                     'derfl' : m(i) = sqrt[ alpha + ---  sum --------------------------- ]
%...                                                     ne  j=1      (z(i+1)-z(i))^2
%...

     choice = 'der1';
%...
%... %%%%%%%%%%%%%%%%%%%%%%%%%%
%... % adapt the initial grid %
%... %%%%%%%%%%%%%%%%%%%%%%%%%%
%...
%... compute the spatial smoothing matrix
%...
     Ainit=diag((1+2*mu)*ones(n+1,1),0)+diag(-mu*ones(n,1),1)+diag(-mu*ones(n,1),-1);
     Ainit(1,1)=Ainit(1,1)-mu;
     Ainit(n+1,n+1)=Ainit(n+1,n+1)-mu;
%...
%... iterative initial grid adaptation
%...
     for iad=1:6
%...
%... compute the monitoring function mon(i)
%...
     mon = monitor(u,z,choice,t);
%...
%... compute the new grid
%...
         suminv=0;
         for i=1:n+1,
             suminv=suminv+1/mon(i);
         end
         for i=1:n+1,
             g(i)=(zR-zL)/(mon(i)*suminv);
         end
         delz=Ainit\g';
         for i=2:n+1,
             z(i)=z(i-1)+delz(i);
         end
%...
%... new initial dependent variables u(n+2,ne)
%...
         for j=1:ne,
             for i=1:n+2,
                 u(i,j) = burgers_exact(z(i),0);
             end
         end         
     end
%...
%... global initial condition
%...
     for j=1:ne,
         x(j:ne+1:n*(ne+1),1)=u(2:n+1,j);
     end
     x(ne+1:ne+1:n*(ne+1),1)=z(2:n+1);
%...
%... %%%%%%%%%%%%%%%%%%%%%%
%... % call to ODE solver %
%... %%%%%%%%%%%%%%%%%%%%%%
%...
     options = odeset('RelTol',1e-6,'AbsTol',1e-6,'Mass',@mass,...
         'MStateDependence','strong','JPattern',JPat,'MvPattern',MvPat,'MassSingular','no','stats','on');
%...     
     [tout, yout] = ode15s(@right_vect,t,x,options);
%...
%... %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%... % visualisation of the solution %
%... %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%...
%... separate the dependent variables and the node positions
%...
     figure(1);
     h1=axes('position',[0.1 0.22 0.8 0.75]);
     h2=axes('position',[0.1 0.1 0.8 0.10]);

     absc=zeros(n+2,1);
     ordo=zeros(n+2,ne);
     ordoth=zeros(n+2,ne);
%...
     axes(h1);
%...
     for i=1:length(tout),
         for j=1:n,
             absc(j+1)=yout(i,(ne+1)*j);
             for k=1:ne,
                 ordo(j+1,k)=yout(i,(ne+1)*j-(ne+1-k));
                 ordoth(j+1,k)=burgers_exact(absc(j+1),tout(i));
             end
         end
         absc(1)=zL;
         absc(n+2)=zR;
         ordo(1,1)=burgers_exact(absc(1),tout(i));
         ordo(n+2,1)=burgers_exact(absc(n+2),tout(i));
         ordoth(1,1)=burgers_exact(absc(1),tout(i));
         ordoth(n+2,1)=burgers_exact(absc(n+2),tout(i));
%...
%...  plot the solution
%...
         plot(absc,ordo(:,1),'k')
         hold on
         plot(absc,ordo(:,1),'.k')
         plot(absc,ordoth(:,1),'r')
     end

     xlabel('z');
     ylabel('u(t,z)');
     axis([0 1 0 1.4]);
%...
     axes(h2);
     for i=1:length(tout),
         for j=1:n,
             absc(j+1)=yout(i,(ne+1)*j);
         end
         absc(1)=zL;
         absc(n+2)=zR;
         plot(absc,((i-1)*1.)*ones(1,n+2),'.k')
         hold on
     end

     xlabel('z');
     ylabel('t');
     tcpu=toc
