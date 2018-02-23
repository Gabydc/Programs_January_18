%... The MatMol Group (2016)
%... Stability analysis of the advection-diffusion-reaction 
%... equation using an explicit Euler IVP solver
%... Dz is kept constant and Dt changes
    close all
    clear all

%... General parameters and spatial grid
    dz  = 0.001;
    z   = 0 : dz : 1;
    n   = length(z);
    nu  = 1;
    k   = -5;
    D   = 0.005;
    dtm = 1/(nu/dz+2*D/(dz^2)-k/2);
    th  = 100*(z-.25);
    th1 = 0:.001:2*pi;

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
    x1 = zeros(1,n);
    for i = 2 : i1
        x1(1)     = 0;
        x1(2:n-1) = (sig+beta)*x0(1:n-2)+(1-sig-2*beta+k*dt)*x0(2:n-1)+beta*x0(3:n);
        x1(n)     = 0;
        x0        = x1;
    end
    x1_Dtmax = x1;

%... At time 0.5
    i2 = find(abs(time-.50) <= dt/2);
    x2 = zeros(1,n);
    for i = i1+1 : i2
        x2(1)     = 0;
        x2(2:n-1) = (sig+beta)*x0(1:n-2)+(1-sig-2*beta+k*dt)*x0(2:n-1)+beta*x0(3:n);
        x2(n)     = 0;
        x0        = x2;
    end
    x2_Dtmax = x2;

%... Eigenvalue loci
    vp_Dtmax = k*dt+(sig+2*beta)*(cos(th1)-1)-j*sig*sin(th1);
    revp_Dtmax = real(vp_Dtmax);
    imvp_Dtmax = imag(vp_Dtmax);
    
%... Plot the results
    figure(1)    
    plot(z,[x0_Dtmax;x1_Dtmax;x2_Dtmax])

%... Clear some workspace variables    
    clear x0 x1 x2 beta sig time dt i1 i2 x0_Dtmax x1_Dtmax x2_Dtmax vp_Dtmax revp_Dtmax imvp_Dtmax



%... SOLUTION WITH Dt = 0.1*Dt_max
%... Initial conditions
    x0       = (0.57)./(exp(th)+exp(-th));
    x0_01Dtmax = x0;

%... Time parameters
    dt   = 0.1*dtm;
    time = 0 : dt : 1;
    sig  = nu*dt/dz;
    beta = D*dt/(dz^2);

%... At time 0.25
    i1 = find(abs(time-.25) <= dt/2);
    x1 = zeros(1,n);
    for i = 2 : i1
        x1(1)     = 0;
        x1(2:n-1) = (sig+beta)*x0(1:n-2)+(1-sig-2*beta+k*dt)*x0(2:n-1)+beta*x0(3:n);
        x1(n)     = 0;
        x0        = x1;
    end
    x1_01Dtmax = x1;

%... At time 0.5
    i2 = find(abs(time-.50) <= dt/2);
    x2 = zeros(1,n);
    for i = i1+1 : i2
        x2(1)     = 0;
        x2(2:n-1) = (sig+beta)*x0(1:n-2)+(1-sig-2*beta+k*dt)*x0(2:n-1)+beta*x0(3:n);
        x2(n)     = 0;
        x0        = x2;
    end
    x2_01Dtmax = x2;

%... Eigenvalue loci
    vp_01Dtmax = k*dt+(sig+2*beta)*(cos(th1)-1)-j*sig*sin(th1);
    revp_01Dtmax = real(vp_01Dtmax);
    imvp_01Dtmax = imag(vp_01Dtmax);

%... Plot the results
    figure(2)    
    plot(z,[x0_01Dtmax;x1_01Dtmax;x2_01Dtmax])

%... Clear some workspace variables
    clear x0 x1 x2 beta sig time dt i1 i2 x0_01Dtmax x1_01Dtmax x2_01Dtmax vp_01Dtmax revp_01Dtmax imvp_01Dtmax


%... SOLUTION WITH Dt = 1.0027*Dt_max
%... Initial conditions
    x0         = (0.57)./(exp(th)+exp(-th));
    x0_10027Dtmax = x0;

%... Time parameters
    dt   = 1.0027*dtm;
    time = 0 : dt : 1;
    sig  = nu*dt/dz;
    beta = D*dt/(dz^2);

%... At time 0.25
    i1 = find(abs(time-.25) <= dt/2);
    x1 = zeros(1,n);
    for i = 2 : i1
        x1(1)     = 0;
        x1(2:n-1) = (sig+beta)*x0(1:n-2)+(1-sig-2*beta+k*dt)*x0(2:n-1)+beta*x0(3:n);
        x1(n)     = 0;
        x0        = x1;
    end
    x1_10027Dtmax = x1;

%... At time 0.5
    i2 = find(abs(time-.50) <= dt/2);
    x2 = zeros(1,n);
    for i = i1+1 : i2
        x2(1)     = 0;
        x2(2:n-1) = (sig+beta)*x0(1:n-2)+(1-sig-2*beta+k*dt)*x0(2:n-1)+beta*x0(3:n);
        x2(n)     = 0;
        x0        = x2;
    end
    x2_10027Dtmax = x2;

%... Eigenvalue loci
    vp_10027Dtmax = k*dt+(sig+2*beta)*(cos(th1)-1)-j*sig*sin(th1);
    revp_10027Dtmax = real(vp_10027Dtmax);
    imvp_10027Dtmax = imag(vp_10027Dtmax);
    
%... Plot the results
    figure(3)
    plot(z,[x0_10027Dtmax;x1_10027Dtmax;x2_10027Dtmax])

%... Clear some workspace variables
    clear x0 x1 x2 beta sig time dt dtm i1 i2 x0_10027Dtmax x1_10027Dtmax x2_10027Dtmax vp_10027Dtmax revp_10027Dtmax imvp_10027Dtmax
