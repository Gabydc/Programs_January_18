%... The MatMol Group (2016)
%... Stability analysis of the advection-diffusion equation using the
%... implicit Euler IVP solver and a downwind finite difference scheme
%... Dt is kept constant and Dz changes
    close all
    clear all

%... General parameters
    nu  = 1;
    D   = 0.005;
    sig = 0.5;
    th1 = 0 : 0.001 : 2*pi;

%... SOLUTION WITH dz = 2*D/(3*nu)
%... Spatial grid
    dz    = 2*D/(3*nu);
    z_2D  = 0 : dz : 1;
    n     = length(z_2D);
    th_2D = 100*(z_2D-0.25);

%... Initial conditions
    x0    = (0.57)./(exp(th_2D)+exp(-th_2D));
    x0_2D = x0;

%... Time parameters
    dt   = dz/(2*nu);
    time = 0 : dt : 1;
    beta = D*dt/(dz^2);

%... At time 0.25
    i1 = find(abs(time-.25) <= dt/2);
    for i = 2 : i1
        x1(1)     = 0;
        x1(2:n-1) = beta*x0(1:n-2)+(1+sig-2*beta)*x0(2:n-1)+(beta-sig)*x0(3:n);
        x1(n)     = 0;
        x0        = x1;
    end
    x1_2D = x1;

%... At time 0.5
    i2 = find(abs(time-.50) <= dt/2);
    for i = i1+1 : i2
        x2(1)     = 0;
        x2(2:n-1) = beta*x0(1:n-2)+(1+sig-2*beta)*x0(2:n-1)+(beta-sig)*x0(3:n);
        x2(n)     = 0;
        x0        = x2;
    end
    x2_2D = x2;

%... Eigenvalue loci
    vp_2D    = (-sig+2*beta)*(cos(th1)-1)-j*sig*sin(th1);
    re_vp_2D = real(vp_2D);
    im_vp_2D = imag(vp_2D);
    
%... Plot the results
    figure(1)
    plot(z_2D,[x0_2D;x1_2D;x2_2D])
    
%... Clear some workspace variables
    clear x0 x1 x2 dz z_2D n th_2D x0_2D x1_2D x2_2D vp_2D re_vp_2D im_vp_2D i i1 i2 dt beta time


%...========================================================
%... SOLUTION WITH dz = 3*D/(3*nu)
%... Spatial grid
    dz    = 3*D/(3*nu);
    z_3D  = 0 : dz : 1;
    n     = length(z_3D);
    th_3D = 100*(z_3D-0.25);

%... Initial conditions
    x0    = (0.57)./(exp(th_3D)+exp(-th_3D));
    x0_3D = x0;

%... Time parameters
    dt   = dz/(2*nu);
    time = 0 : dt : 1;
    beta = D*dt/(dz^2);

%... At time 0.25
    i1 = find(abs(time-.25) <= dt/2);
    for i = 2 : i1
        x1(1)     = 0;
        x1(2:n-1) = beta*x0(1:n-2)+(1+sig-2*beta)*x0(2:n-1)+(beta-sig)*x0(3:n);
        x1(n)     = 0;
        x0        = x1;
    end
    x1_3D = x1;

%... At time 0.5
    i2 = find(abs(time-.50) <= dt/2);
    for i = i1+1 : i2
        x2(1)     = 0;
        x2(2:n-1) = beta*x0(1:n-2)+(1+sig-2*beta)*x0(2:n-1)+(beta-sig)*x0(3:n);
        x2(n)     = 0;
        x0        = x2;
    end
    x2_3D = x2;

%... Eigenvalue loci
    vp_3D    = (-sig+2*beta)*(cos(th1)-1)-j*sig*sin(th1);
    re_vp_3D = real(vp_3D);
    im_vp_3D = imag(vp_3D);
%... Plot the results
    figure(2)
    plot(z_3D,[x0_3D;x1_3D;x2_3D])
%... Clear some workspace variables
    clear x0 x1 x2 dz z_3D n th_3D x0_3D x1_3D x2_3D vp_3D re_vp_3D im_vp_3D i i1 i2 dt beta time




%... SOLUTION WITH dz = 4*D/(3*nu)
%... Spatial grid
    dz    = 4*D/(3*nu);
    z_4D  = 0 : dz : 1;
    n     = length(z_4D);
    th_4D = 100*(z_4D-0.25);

%... Initial conditions
    x0    = (0.57)./(exp(th_4D)+exp(-th_4D));
    x0_4D = x0;

%... Time parameters
    dt   = dz/(2*nu);
    time = 0 : dt : 1;
    beta = D*dt/(dz^2);

%... At time 0.25
    i1 = find(abs(time-.25) <= dt/2);
    for i = 2 : i1
        x1(1)     = 0;
        x1(2:n-1) = beta*x0(1:n-2)+(1+sig-2*beta)*x0(2:n-1)+(beta-sig)*x0(3:n);
        x1(n)     = 0;
        x0        = x1;
    end
    x1_4D = x1;

%... At time 0.5
    i2 = find(abs(time-.50) <= dt/2);
    for i = i1+1 : i2
        x2(1)     = 0;
        x2(2:n-1) = beta*x0(1:n-2)+(1+sig-2*beta)*x0(2:n-1)+(beta-sig)*x0(3:n);
        x2(n)     = 0;
        x0        = x2;
    end
    x2_4D = x2;

%... Eigenvalue loci
    vp_4D    = (-sig+2*beta)*(cos(th1)-1)-j*sig*sin(th1);
    re_vp_4D = real(vp_4D);
    im_vp_4D = imag(vp_4D);
%... Plot the results
    figure(3)
    plot(z_4D,[x0_4D;x1_4D;x2_4D])
%... Clear some workspace variables
    clear x0 x1 x2 dz z_4D n th_4D x0_4D x1_4D x2_4D vp_4D re_vp_4D im_vp_4D i i1 i2 dt beta time

