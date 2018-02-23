%... The MatMol Group (2016)
%... Solution of the advection equation using the Euler IVP
%... solver and a finite difference scheme for different values of
%... sigma
    close all
    clear all

%... Spatial grid
    dz  = 0.01;
    z   = 0 : dz : 1;
    n   = length(z);
    nu  = 1;
    th  = 100*(z-.25);
    th1 = 0 : 0.001 : 2*pi;

%... SOLUTION WITH SIGMA = 0.95
%... Initial conditions
    x0_s095 = (.57)./(exp(th)+exp(-th));

%... Time parameters
    sig  = 0.95;
    dt   = sig*dz/nu;
    time = 0 : dt : 1;

%... At time 0.25
    i1 = find(abs(time-0.25) <= dt/2);
    xk = x0_s095;
    for ii = 2 : i1
        x1_s095(1)   = 0;
        x1_s095(2:n) = (1-sig)*xk(2:n)+sig*xk(1:n-1);
        xk           = x1_s095;
    end

%... At time 0.5
    i2   = find(abs(time-0.50) <= dt/2);
    for ii = i1+1:i2
        x2_s095(1)   = 0;
        x2_s095(2:n) = (1-sig)*xk(2:n)+sig*xk(1:n-1);
        xk           = x2_s095;
    end
    figure
    x_dz001_s095 = [x0_s095;x1_s095;x2_s095];
    plot(z,x_dz001_s095)

%... Eigenvalue loci
    vp_s095   = sig*(cos(th1)-1-j*sin(th1));
    revp_s095 = real(vp_s095);
    imvp_s095 = imag(vp_s095);
    
%... SOLUTION WITH SIGMA = 1.00
%... Initial conditions
    x0_s100 = (.57)./(exp(th)+exp(-th));

%... Time parameters
    sig  = 1.00;
    dt   = sig*dz/nu;
    time = 0 : dt : 1;

%... At time 0.25
    i1 = find(abs(time-0.25) <= dt/2);
    xk = x0_s100;
    for ii = 2 : i1
        x1_s100(1)   = 0;
        x1_s100(2:n) = (1-sig)*xk(2:n)+sig*xk(1:n-1);
        xk           = x1_s100;
    end

%... At time 0.5
    i2   = find(abs(time-0.50) <= dt/2);
    for ii = i1+1:i2
        x2_s100(1)   = 0;
        x2_s100(2:n) = (1-sig)*xk(2:n)+sig*xk(1:n-1);
        xk           = x2_s100;
    end
    figure
    x_dz001_s100 = [x0_s100;x1_s100;x2_s100];
    plot(z,x_dz001_s100)

%... Eigenvalue loci
    vp_s100   = sig*(cos(th1)-1-j*sin(th1));
    revp_s100 = real(vp_s100);
    imvp_s100 = imag(vp_s100);


%... SOLUTION WITH SIGMA = 1.05
%... Initial conditions
    x0_s105 = (.57)./(exp(th)+exp(-th));

%... Time parameters
    sig  = 1.05;
    dt   = sig*dz/nu;
    time = 0 : dt : 1;

%... At time 0.25
    i1 = find(abs(time-0.25) <= dt/2);
    xk = x0_s105;
    for ii = 2 : i1
        x1_s105(1)   = 0;
        x1_s105(2:n) = (1-sig)*xk(2:n)+sig*xk(1:n-1);
        xk           = x1_s105;
    end

%... At time 0.5
    i2   = find(abs(time-0.50) <= dt/2);
    for ii = i1+1:i2
        x2_s105(1)   = 0;
        x2_s105(2:n) = (1-sig)*xk(2:n)+sig*xk(1:n-1);
        xk           = x2_s105;
    end
    figure
    x_dz001_s105 = [x0_s105;x1_s105;x2_s105];
    plot(z,x_dz001_s105)

%... Eigenvalue loci
    vp_s105   = sig*(cos(th1)-1-j*sin(th1));
    revp_s105 = real(vp_s105);
    imvp_s105 = imag(vp_s105);
