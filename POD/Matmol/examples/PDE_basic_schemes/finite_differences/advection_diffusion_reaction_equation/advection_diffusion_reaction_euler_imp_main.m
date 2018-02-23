%... The MatMol Group (2016)
%... Stability analysis of the advection-diffusion-reaction 
%... equation using an implicit Euler IVP solver
%... Dz is kept constant and Dt changes
    close all
    clear all

%... General parameters and spatial grid
    dz        = 0.001;
    z         = 0:dz:1;
    n         = length(z);
    nu        = 1;
    D         = 0.005;
    dtm       = 1/(nu/dz+2*D/(dz^2));
    th(1:n,1) = 100*(z-.25);
    th1       = 0:.001:2*pi;

%... SOLUTION WITH Dt = Dt_max
%... Initial conditions
    x0       = (0.57)./(exp(th)+exp(-th));
    x0_Dtmax = x0;

%... Time parameters
    dt   = dtm;
    time = 0 : dt : 1;
    sig  = nu*dt/dz;
    beta = D*dt/(dz^2);

%... At time 0.25
    i1 = find(abs(time-.25) <= dt/2);
    x1 = zeros(n,1);
    M = diag((1+sig+2*beta)*ones(1,n-2),0)+diag(-beta*ones(1,n-3),1)+diag(-(sig+beta)*ones(1,n-3),-1);
    for i = 1 : i1
        x1(1,1)     = 0;
        x1(2:n-1,1) = linsolve(M,x0(2:n-1,1));
        x1(n,1)     = 0;
        x0          = x1;
    end
    x1_Dtmax = x1;

%... At time 0.5
    i2 = find(abs(time-.50) <= dt/2);
    x2 = zeros(n,1);
    for i = i1+1 : i2
        x2(1,1)     = 0;
        x2(2:n-1,1) = linsolve(M,x0(2:n-1,1));
        x2(n,1)     = 0;
        x0          = x2;
    end
    x2_Dtmax = x2;

%... Plot the results
    figure(1)
    plot(z,[x0_Dtmax x1_Dtmax x2_Dtmax])    

%... Clear some workspace variables
    clear x0 x1 x2 i1 i2 dt time sig beta x0_Dtmax x1_Dtmax x2_Dtmax


%... SOLUTION WITH Dt = 10*Dt_max
%... Initial conditions
    x0       = (0.57)./(exp(th)+exp(-th));
    x0_10Dtmax = x0;

%... Time parameters
    dt   = 10*dtm;
    time = 0 : dt : 1;
    sig  = nu*dt/dz;
    beta = D*dt/(dz^2);

%... At time 0.25
    i1 = find(abs(time-.25) <= dt/2);
    x1 = zeros(n,1);
    M = diag((1+sig+2*beta)*ones(1,n-2),0)+diag(-beta*ones(1,n-3),1)+diag(-(sig+beta)*ones(1,n-3),-1);
    for i = 1 : i1
        x1(1,1)     = 0;
        x1(2:n-1,1) = linsolve(M,x0(2:n-1,1));
        x1(n,1)     = 0;
        x0          = x1;
    end
    x1_10Dtmax = x1;

%... At time 0.5
    i2 = find(abs(time-.50) <= dt/2);
    x2 = zeros(n,1);
    for i = i1+1 : i2
        x2(1,1)     = 0;
        x2(2:n-1,1) = linsolve(M,x0(2:n-1,1));
        x2(n,1)     = 0;
        x0          = x2;
    end
    x2_10Dtmax = x2;
    
%... Plot the results
    figure(2)
    plot(z,[x0_10Dtmax x1_10Dtmax x2_10Dtmax])    

%... Clear some workspace variables
    clear x0 x1 x2 i1 i2 dt time sig beta x0_10Dtmax x1_10Dtmax x2_10Dtmax



%... SOLUTION WITH Dt = 100*Dt_max
%... Initial conditions
    x0          = (0.57)./(exp(th)+exp(-th));
    x0_100Dtmax = x0;

%... Time parameters
    dt   = 100*dtm;
    time = 0 : dt : 1;
    sig  = nu*dt/dz;
    beta = D*dt/(dz^2);

%... At time 0.25
    i1 = find(abs(time-.25) <= dt/2);
    x1 = zeros(n,1);
    M = diag((1+sig+2*beta)*ones(1,n-2),0)+diag(-beta*ones(1,n-3),1)+diag(-(sig+beta)*ones(1,n-3),-1);
    for i = 1 : i1
        x1(1,1)     = 0;
        x1(2:n-1,1) = linsolve(M,x0(2:n-1,1));
        x1(n,1)     = 0;
        x0          = x1;
    end
    x1_100Dtmax = x1;

%... At time 0.5
    i2 = find(abs(time-.50) <= dt/2);
    x2 = zeros(n,1);
    for i = i1+1 : i2
        x2(1,1)     = 0;
        x2(2:n-1,1) = linsolve(M,x0(2:n-1,1));
        x2(n,1)     = 0;
        x0          = x2;
    end
    x2_100Dtmax = x2;

%... Plot the results
    figure(3)
    plot(z,[x0_100Dtmax x1_100Dtmax x2_100Dtmax])

%... Clear some workspace variables
    clear x0 x1 x2 i1 i2 dt time sig beta x0_100Dtmax x1_100Dtmax x2_100Dtmax

