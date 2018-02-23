%... The MatMol Group (2016)
%... Solution of the advection equation using the Leap frog IVP
%... solver and a finite difference scheme for different values of
%... sigma
    close all
    clear all

%... Spatial grid
    dz  = 0.01;
    z1  = 0 : dz : 1;
    n   = length(z1);
    nu  = 1;
    th1 = 100*(z1-.25);

%... SOLUTION WITH Dz = 0.01 and sigma = 1.00
%... Initial conditions
    x0      = (.57)./(exp(th1)+exp(-th1));
    x0_s100 = x0;
%... Time parameters 
    sig  = 1;
    dt   = sig*dz/nu;
    time = 0 : dt : 1;

%... At time 0.25
    i1 = find(abs(time-.25) <= dt/2);
    x1(1)   = 0;
    x1(2:n) = (1-sig)*x0(2:n)+sig*x0(1:n-1);
    x2      = zeros(1,n);
    for i = 3 : i1
        x2(1)     = 0;
        x2(2:n-1) = x0(2:n-1)-sig*(x1(3:n)-x1(1:n-2));
        x2(n)     = 0;
        x0        = x1;
        x1        = x2;
    end
    x2_s100 = x2;

%... At time 0.5
    i2 = find(abs(time-.50) <= dt/2);
    x3 = zeros(1,n);
    for i = i1+1 : i2
        x3(1)     = 0;
        x3(2:n-1) = x0(2:n-1)-sig*(x1(3:n)-x1(1:n-2));
        x3(n)     = 0;
        x0        = x1;
        x1        = x3;
    end
    x3_s100 = x3;
    clear x0 x1 x2 x3 sig dt time i1 i2
    figure
    x_dz001_s100 = [x0_s100;x2_s100;x3_s100];
    plot(z1,x_dz001_s100)


%... SOLUTION WITH Dz = 0.01 and sigma = 0.5
%... Initial conditions
    x0      = (.57)./(exp(th1)+exp(-th1));
    x0_s050 = x0;
%... Time parameters 
    sig  = 0.5;
    dt   = sig*dz/nu;
    time = 0 : dt : 1;

%... At time 0.25
    i1 = find(abs(time-.25) <= dt/2);
    x1(1)   = 0;
    x1(2:n) = (1-sig)*x0(2:n)+sig*x0(1:n-1);
    x2      = zeros(1,n);
    for i = 3 : i1
        x2(1)     = 0;
        x2(2:n-1) = x0(2:n-1)-sig*(x1(3:n)-x1(1:n-2));
        x2(n)     = 0;
        x0        = x1;
        x1        = x2;
    end
    x2_s050 = x2;

%... At time 0.5
    i2 = find(abs(time-.50) <= dt/2);
    x3 = zeros(1,n);
    for i = i1+1 : i2
        x3(1)     = 0;
        x3(2:n-1) = x0(2:n-1)-sig*(x1(3:n)-x1(1:n-2));
        x3(n)     = 0;
        x0        = x1;
        x1        = x3;
    end
    x3_s050 = x3;
    clear x0 x1 x2 x3 sig dt time i1 i2
    figure
    x_dz001_s050 = [x0_s050;x2_s050;x3_s050];
    plot(z1,x_dz001_s050)


%... SOLUTION WITH Dz = 0.01 and sigma = 1.003
%... Initial conditions
    x0      = (.57)./(exp(th1)+exp(-th1));
    x0_s1003 = x0;
%... Time parameters 
    sig  = 1.003;
    dt   = sig*dz/nu;
    time = 0 : dt : 1;

%... At time 0.25
    i1 = find(abs(time-.25) <= dt/2);
    x1(1)   = 0;
    x1(2:n) = (1-sig)*x0(2:n)+sig*x0(1:n-1);
    x2      = zeros(1,n);
    for i = 3 : i1
        x2(1)     = 0;
        x2(2:n-1) = x0(2:n-1)-sig*(x1(3:n)-x1(1:n-2));
        x2(n)     = 0;
        x0        = x1;
        x1        = x2;
    end
    x2_s1003 = x2;

%... At time 0.5
    i2 = find(abs(time-.50) <= dt/2);
    x3 = zeros(1,n);
    for i = i1+1 : i2
        x3(1)     = 0;
        x3(2:n-1) = x0(2:n-1)-sig*(x1(3:n)-x1(1:n-2));
        x3(n)     = 0;
        x0        = x1;
        x1        = x3;
    end
    x3_s1003 = x3;
    clear x0 x1 x2 x3 sig dt time i1 i2
    figure
    x_dz001_s1003 = [x0_s1003;x2_s1003;x3_s1003];
    plot(z1,x_dz001_s1003)

    
%... Spatial grid
    dz  = 0.001;
    z2  = 0 : dz : 1;
    n   = length(z2);
    nu  = 1;
    th2 = 100*(z2-.25);

%... SOLUTION WITH Dz = 0.001 and sigma = 1.00
%... Initial conditions
    x0      = (.57)./(exp(th2)+exp(-th2));
    x0_s100 = x0;
%... Time parameters 
    sig  = 1;
    dt   = sig*dz/nu;
    time = 0 : dt : 1;

%... At time 0.25
    i1 = find(abs(time-.25) <= dt/2);
    x1(1)   = 0;
    x1(2:n) = (1-sig)*x0(2:n)+sig*x0(1:n-1);
    x2      = zeros(1,n);
    for i = 3 : i1
        x2(1)     = 0;
        x2(2:n-1) = x0(2:n-1)-sig*(x1(3:n)-x1(1:n-2));
        x2(n)     = 0;
        x0        = x1;
        x1        = x2;
    end
    x2_s100 = x2;

%... At time 0.5
    i2 = find(abs(time-.50) <= dt/2);
    x3 = zeros(1,n);
    for i = i1+1 : i2
        x3(1)     = 0;
        x3(2:n-1) = x0(2:n-1)-sig*(x1(3:n)-x1(1:n-2));
        x3(n)     = 0;
        x0        = x1;
        x1        = x3;
    end
    x3_s100 = x3;
    clear x0 x1 x2 x3 sig dt time i1 i2
    figure
    x_dz0001_s100 = [x0_s100;x2_s100;x3_s100];
    plot(z2,x_dz0001_s100)


%... SOLUTION WITH Dz = 0.001 and sigma = 0.5
%... Initial conditions
    x0      = (.57)./(exp(th2)+exp(-th2));
    x0_s050 = x0;
%... Time parameters 
    sig  = 0.5;
    dt   = sig*dz/nu;
    time = 0 : dt : 1;

%... At time 0.25
    i1 = find(abs(time-.25) <= dt/2);
    x1(1)   = 0;
    x1(2:n) = (1-sig)*x0(2:n)+sig*x0(1:n-1);
    x2      = zeros(1,n);
    for i = 3 : i1
        x2(1)     = 0;
        x2(2:n-1) = x0(2:n-1)-sig*(x1(3:n)-x1(1:n-2));
        x2(n)     = 0;
        x0        = x1;
        x1        = x2;
    end
    x2_s050 = x2;

%... At time 0.5
        i2 = find(abs(time-.50) <= dt/2);
    x3 = zeros(1,n);
    for i = i1+1 : i2
        x3(1)     = 0;
        x3(2:n-1) = x0(2:n-1)-sig*(x1(3:n)-x1(1:n-2));
        x3(n)     = 0;
        x0        = x1;
        x1        = x3;
    end
    x3_s050 = x3;
    clear x0 x1 x2 x3 sig dt time i1 i2
    figure
    x_dz0001_s050 = [x0_s050;x2_s050;x3_s050];
    plot(z2,x_dz0001_s050)


%... SOLUTION WITH Dz = 0.001 and sigma = 1.0012
%... Initial conditions
    x0        = (.57)./(exp(th2)+exp(-th2));
    x0_s10012 = x0;
%... Time parameters 
    sig  = 1.0012;
    dt   = sig*dz/nu;
    time = 0 : dt : 1;

%... At time 0.25
    i1 = find(abs(time-.25) <= dt/2);
    x1(1)   = 0;
    x1(2:n) = (1-sig)*x0(2:n)+sig*x0(1:n-1);
    x2      = zeros(1,n);
    for i = 3 : i1
        x2(1)     = 0;
        x2(2:n-1) = x0(2:n-1)-sig*(x1(3:n)-x1(1:n-2));
        x2(n)     = 0;
        x0        = x1;
        x1        = x2;
    end
    x2_s10012 = x2;

%... At time 0.5
    i2 = find(abs(time-.50) <= dt/2);
    x3 = zeros(1,n);
    for i = i1+1 : i2
        x3(1)     = 0;
        x3(2:n-1) = x0(2:n-1)-sig*(x1(3:n)-x1(1:n-2));
        x3(n)     = 0;
        x0        = x1;
        x1        = x3;
    end
    x3_s10012 = x3;
    clear x0 x1 x2 x3 sig dt time i1 i2
    figure
    x_dz0001_s10012 = [x0_s10012;x2_s10012;x3_s10012];
    plot(z2,x_dz0001_s10012)



