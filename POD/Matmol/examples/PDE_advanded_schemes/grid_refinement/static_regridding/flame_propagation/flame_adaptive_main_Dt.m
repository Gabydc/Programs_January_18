%... The MatMol Group (2016)
%...
%...  Flame propagation problem
%...
%... The classical Dwyer's and Sanders' flame propagation problem:
%...
%...    rho  = rho   - NDA*rho                                 (1)
%...       t      zz    
%...
%...    T  = T   + NDA*rho                                     (2)
%...     t    zz    
%...
%... where NDA = 3.52e6*exp(-4/T) and 
%...
%... The boundary conditions are given by
%...
%... rho (0,t) = 0                rho (1,t) = 0                (3)
%...    z                            z
%...
%... T (0,t) = 0                  T(1,t) = f(t)                (4)
%...  z
%...
%... where f(t) = 0.2 + t/2e-4  for t < 2e-4 and f(t) = 1.2 otherwise
%...
%... The initial conditions are given by
%...
%... rho(z,0) = 1                 T(z,0) = 0.2                  (5)
%...
%... The following code computes a solution to eqs. (1-5)
%...
     close all
     clear all
%...
%... Start a stopwatch timer
     tic
%...
%... Set global variables
     global z0 zL nz D1
     global nsteps maxsteps
%...
%... Spatial grid
     z0 = 0.0;
     zL = 1.0;
     nz = 201;
     dz = (zL-z0)/(nz-1);
     z = [z0:dz:zL]';
%...
%... Initial conditions
     rho = ones(nz,1);
     T = 0.2*ones(nz,1);
     x(1:2:2*nz-1) = rho;
     x(2:2:2*nz) = T;
     x = x';
%...
%... parameters of the adaptive grid
     npdes=2;
     nzmax=1001;
     alpha=0.05;
     beta=100;
     tolz=0.02;
     bound=1.4;
     imesh=1;
     ilim=0;
%...
%... differentiation matrix
     D1 = three_point_centered_D1(z);
%...  
%... call to ODE solver 
%...     
     t0 = 0;
     tf = 0.006;
     yout = x;
     zout = z;
     nzout = nz;
     tout = t0;
%...
%... solver to stop after this many steps:
     maxsteps = 10;
%...
%... initial situation
%...
     figure(1)
     subplot('position',[0.1 0.3 0.8 0.6])
     plot(z,rho,'b');
%     xlabel('z')
     ylabel('\rho(z,t)');
%     title('Flame propagation')
     axis([0 1 0 1.5])
     hold on
     subplot('position',[0.1 0.08 0.8 0.17])
     plot(z,1000*t0*ones(nz,1),'.k')
     xlabel('z');
     ylabel('t')
     axis([0 1 0 6])
     hold on
%...     
     figure(2)
     subplot('position',[0.1 0.3 0.8 0.6])
     plot(z,T,'r');
%     xlabel('z')
     ylabel('T(z,t)');
%     title('Flame propagation')
     axis([0 1 0 1.5])
     hold on
     subplot('position',[0.1 0.08 0.8 0.17])
     plot(z,1000*t0*ones(nz,1),'.k')
     xlabel('z');
     ylabel('t')
     axis([0 1 0 6])
     hold on
%...
%... capture first frame
%     iframe = 1;
%     mov(iframe) = getframe(1);
%...
%...     
     Dt = tf/10;
     tk = t0;
     tkp1 = tk + Dt;
     tspan = [tk tkp1]
%...     
cont = 1;
zz{cont}   = z;
TT{cont}   = T';
rhor{cont} = rho';
tt(cont)   = t0;
     while tkp1 <= tf
         cont = cont + 1;
%...
%... do the integration for maxsteps steps in a loop until t becomes larger than tf
        options = odeset('RelTol',1e-3,'AbsTol',1e-3);
        options = odeset(options,'JPattern',JPattern(nz));
        [t,y] = ode23s(@flame_adaptive_pdes,tspan,x,options);
%...       
        tk = t(end);
        tkp1 = tk + Dt;
        tspan = [tk tkp1];
        x = [];
        x = y(end,:);
%...
%... refine the grid
%...
        [z_new,nz_new,ier,tolz_new] = agereg(z,x,npdes,nzmax,alpha,beta,tolz,bound,imesh,ilim);
%...
%... interpolate the dependent variables
%...
        rho = x(1:2:2*nz-1);
        T = x(2:2:2*nz);
        rho_new = spline(z,rho,z_new);
        T_new = spline(z,T,z_new);
        x = [];
        x(1:2:2*nz_new-1) = rho_new;
        x(2:2:2*nz_new) = T_new;
        x = x';
        yout = [yout ; x];
%...
        z = z_new';
        nz = nz_new;
        tolz = tolz_new;
        zout = [zout ; z];
        nzout = [nzout ; nz];
        tout = [tout ; tk];
%...
%... plot intermediate results
        figure(1)
        subplot('position',[0.1 0.3 0.8 0.6])
        plot(z,rho_new,'b');
        axis([0 1 0 1.2])
%        hold on
%        plot(z,rho_new,'.b',z,T_new,'.r');
%...
         subplot('position',[0.1 0.08 0.8 0.17])
         plot(z,1000*tk*ones(nz,1),'.k')
%...        
        figure(2)
        subplot('position',[0.1 0.3 0.8 0.6])
        plot(z,T_new,'r');
%        hold on
%        plot(z,rho_new,'.b',z,T_new,'.r');
%...
         subplot('position',[0.1 0.08 0.8 0.17])
         plot(z,1000*tk*ones(nz,1),'.k')
         
         zz{cont}   = z;
         TT{cont}   = T_new;
         rhor{cont} = rho_new;
         tt(cont)   = tk;
%...
%... capture graph in a movie        
%        iframe = iframe+1;
%        mov(iframe) = getframe(1);
%...
%... compute a new differentiation matrix
%...
         D1 = three_point_centered_D1(z);
%...
     end
%...
%... read the stopwatch timer
%...
    nav = sum(nzout)/length(nzout)
    tcpu=toc;
     
