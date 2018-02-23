%... The MatMol Group (2016)
%...
%...  Korteweg-de Vries equation
%...
%... The classical Korteweg-de Vries equation:
%...
%...    u  = -6*u*u  - u                                  (1)
%...     t         z    zzz
%...
%... The propagation of a single soliton can be written as  
%...
%... u(z,t)= 0.5*s*sech(0.5*sqrt(s)*(z-s*t))^2            (2)
%...
%... The initial conditions (ICs) for eq. (1) is derived
%... from (2)
%...
%... The following code computes a solution to eqs. (1-3)
%...
     close all
     clear all
%...
%... Start a stopwatch timer
     tic
%...
%... Set global variables
     global s
     global z0 zL nz D1
     global nsteps maxsteps tprint tflag tout
%...
%... Spatial grid
     z0=-30.0;
     zL=70.0;
     nz=201;
     nzout=nz;
     dz=(zL-z0)/(nz-1);
     z=[z0:dz:zL]';
%...
%... Initial conditions
     s=0.5;
     x=kdv3_exact(z,0);
%...
%... parameters of the adaptive grid
     npdes=1;
     nzmax=1001;
     alpha=0;
     beta=100;
     tolz=0.005;
     bound=1.1;
     imesh=0;
     ilim=0;
%...
%... refine the initial grid
%...
     [z_new,nz_new,ier,tolz_new]=agereg(z,x,npdes,nzmax,alpha,beta,tolz,bound,imesh,ilim);
%...
%... interpolate the dependent variables
%...
     x_new=spline(z,x,z_new);
     x=x_new';
     z=z_new';
     nz=nz_new;
     tolz=tolz_new;
%...
%... differentiation matrix
     D1=three_point_centered_D1(z);
%...  
%... call to ODE solver 
%...     
     t0=0;
     tf=100;
     dt = 5;
     yout=x;
     zout=z;
     nzout=[nzout ; nz];
     tout=t0;
%...
%... solver to stop after this many steps:
     maxsteps = 10;
%...
%... initial situation
%...
     figure(1)
     subplot('position',[0.1 0.3 0.8 0.6])
     plot(z,x);
     ylabel('x(z,t)');
%     title('Korteweg-de Vries equation')
     axis([-30 70 0 0.3])
     hold on
     subplot('position',[0.1 0.08 0.8 0.17])
     plot(z,t0*ones(nz,1),'.b')
     ylabel('t');
     xlabel('z');
     axis([-30 70 0 tf])
     hold on
%...     
     tk = t0;
     tspan = [t0 tf];
     tprint = dt;
     count = 1;
     while tk <= tf-1.e-5
         if count == 1
             xnum{count} = x;
             xana{count} = x;
             zmat{count} = z;
             time(count)   = tk;
         else
             xnum{count} = x;
             xana{count} = yexact;
             zmat{count} = z;
             time(count) = tk;
         end
     count = count + 1;
%...
%... initialize step counter
        nsteps = 0;
%...
%... do the integration for maxsteps steps in a loop until t becomes larger
%... than tf
        options = odeset('RelTol',1e-3,'AbsTol',1e-6);
        options = odeset(options,'Events',@count_steps);
        options = odeset(options,'JPattern',jpattern(nz));
        [t,y,te,ye,ie] = ode23s(@kdv3_adaptive_pde,tspan,x,options);
%...       
        tk=t(end);
        tspan=[tk tf];
        x=[];
        x=y(end,:);
%...
%... refine the grid
%...
        [z_new,nz_new,ier,tolz_new]=agereg(z,x,npdes,nzmax,alpha,beta,tolz,bound,imesh,ilim);
%...
%... interpolate the dependent variables
%...
        x_new=spline(z,x,z_new);
        x=x_new';
        yout=[yout ; x];
%...
        z=z_new';
        nz=nz_new;
        tolz=tolz_new;
        zout=[zout ; z];
        nzout=[nzout ; nz];
        tout=[tout ; tk];
%...
%... plot intermediate results
        if tflag >= 0
            figure(1)
            subplot('position',[0.1 0.3 0.8 0.6])
            plot(z,x);
%...     
            yexact=kdv3_exact(z,tk);
            plot(z,yexact(1:length(z)),'r')
%...
            subplot('position',[0.1 0.08 0.8 0.17])
            plot(z,tk*ones(nz,1),'.b')
            tprint = tprint + dt;
        end
%...
%... compute a new differentiation matrix
%...
         D1=three_point_centered_D1(z);
%...
     end
%...
%... read the stopwatch timer
     tcpu=toc
     nav = sum(nzout)/length(nzout)
     