%... The MatMol Group (2016)
%... Stability analysis of the advection-diffusion equation using the
%... explicit Euler IVP solver and a downwind finite difference scheme
%... Dz is kept constant and Dt changes
    close all
    clear all

%... General parameters and spatial grid
    nu  = 1;
    D   = 0.005;
    dz  = 0.001;
    z   = 0 : dz : 1;
    n   = length(z);
    th  = 100*(z-.25);
    th1 = 0:.001:2*pi;
    dtm = 1e-3/9;

%... SOLUTION WITH Dt = Dt_max
%... Initial conditions
    x0       = (0.57)./(exp(th)+exp(-th));
    x0_Dtmax = x0;

%... Time parameters
    dt   = dtm;
    sig  = nu*dt/dz;
    time = 0:dt:1;
    beta = D*dt/(dz^2);

%... At time 0.25
    i1 = find(abs(time-.25) <= dt/2);
    x1 = zeros(1,n);
    for ii = 2 : i1
        x1(1)     = 0;
        x1(2:n-1) = beta*x0(1:n-2)+(1+sig-2*beta)*x0(2:n-1)+(beta-sig)*x0(3:n);
        x1(n)     = 0;
        x0        = x1;
    end
    x1_Dtmax = x1;

%... At time 0.5
    i2 = find(abs(time-.50) <= dt/2);
    x2 = zeros(1,n);
    for ii = i1+1 : i2
        x2(1)     = 0;
        x2(2:n-1) = beta*x0(1:n-2)+(1+sig-2*beta)*x0(2:n-1)+(beta-sig)*x0(3:n);
        x2(n)     = 0;
        x0        = x2;
    end
    x2_Dtmax = x2;

%... Eigenvalue loci
    vp_Dtmax = (-sig+2*beta)*(cos(th1)-1)-j*sig*sin(th1);
    revp_Dtmax = real(vp_Dtmax);
    imvp_Dtmax = imag(vp_Dtmax);
%... Plot the results    
    figure(1)
    plot(z,[x0_Dtmax;x1_Dtmax;x2_Dtmax])
%... Clear some workspace variables    
    clear x1_Dtmax x2_Dtmax x0 x1 x2 beta dt sig time i1 i2 vp_Dtmax revp_Dtmax imvp_Dtmax ii 


%... SOLUTION WITH Dt = 0.5*Dt_max
%... Initial conditions
    x0 = (0.57)./(exp(th)+exp(-th));

%... Time parameters
    dt   = 0.5*dtm;
    sig  = nu*dt/dz;
    time = 0:dt:1;
    beta = D*dt/(dz^2);

%... At time 0.25
    i1 = find(abs(time-.25) <= dt/2);
    x1 = zeros(1,n);
    for ii = 2 : i1
        x1(1)     = 0;
        x1(2:n-1) = beta*x0(1:n-2)+(1+sig-2*beta)*x0(2:n-1)+(beta-sig)*x0(3:n);
        x1(n)     = 0;
        x0        = x1;
    end
    x1_05Dtmax = x1;

%... At time 0.5
    i2 = find(abs(time-.50) <= dt/2);
    x2 = zeros(1,n);
    for ii = i1+1 : i2
        x2(1)     = 0;
        x2(2:n-1) = beta*x0(1:n-2)+(1+sig-2*beta)*x0(2:n-1)+(beta-sig)*x0(3:n);
        x2(n)     = 0;
        x0        = x2;
    end
    x2_05Dtmax = x2;

%... Eigenvalue loci
    vp_05Dtmax = (-sig+2*beta)*(cos(th1)-1)-j*sig*sin(th1);
    revp_05Dtmax = real(vp_05Dtmax);
    imvp_05Dtmax = imag(vp_05Dtmax);
%... Plot the results    
    figure(2)
    plot(z,[x0_Dtmax;x1_05Dtmax;x2_05Dtmax])
%... Clear some workspace variables    
    clear x1_05Dtmax x2_05Dtmax x0 x1 x2 beta dt sig time i1 i2 vp_05Dtmax revp_05Dtmax imvp_05Dtmax ii 


%... SOLUTION WITH Dt = 1.0033*Dt_max
%... Initial conditions
    x0       = (0.57)./(exp(th)+exp(-th));

%... Time parameters
    dt   = 1.0033*dtm;
    sig  = nu*dt/dz;
    time = 0:dt:1;
    beta = D*dt/(dz^2);

%... At time 0.25
    i1 = find(abs(time-.25) <= dt/2);
    x1 = zeros(1,n);
    for ii = 2 : i1
        x1(1)     = 0;
        x1(2:n-1) = beta*x0(1:n-2)+(1+sig-2*beta)*x0(2:n-1)+(beta-sig)*x0(3:n);
        x1(n)     = 0;
        x0        = x1;
    end
    x1_10033Dtmax = x1;

%... At time 0.5
    i2 = find(abs(time-.50) <= dt/2);
    x2 = zeros(1,n);
    for ii = i1+1 : i2
        x2(1)     = 0;
        x2(2:n-1) = beta*x0(1:n-2)+(1+sig-2*beta)*x0(2:n-1)+(beta-sig)*x0(3:n);
        x2(n)     = 0;
        x0        = x2;
    end
    x2_10033Dtmax = x2;

%... Eigenvalue loci
    vp_10033Dtmax = (-sig+2*beta)*(cos(th1)-1)-j*sig*sin(th1);
    revp_10033Dtmax = real(vp_10033Dtmax);
    imvp_10033Dtmax = imag(vp_10033Dtmax);
%... Plot the results    
    figure(3)
    plot(z,[x0_Dtmax;x1_10033Dtmax;x2_10033Dtmax])
%... Clear some workspace variables    
    clear x1_10033Dtmax x2_10033Dtmax x0 x1 x2 beta dt sig time i1 i2 vp_10033Dtmax revp_10033Dtmax imvp_10033Dtmax ii 

