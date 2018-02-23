%...  The MatMol Group (2016)
%...
%... The well known Buckley-Leverett equation:
%...
%...   xt + (f(x))z = eps*[nu(x)*xz]z
%...
%...   with    f(x) = (x^2)/(x^2+(1-x)^2)
%...
%...           nu(x) = 4x*(1-x)
%...
%...           0 < t < 1
%...
%...   ICs:  x(z,0) = 1-3z   0 < z < 1/3
%...                  0      1/3 < z < 1
%...
%...
%...   BCs : x(0,t) = 1
%...         x(1,t) = 0
%...
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
%...  nu      diffusion coefficient
%...
%...
%... the following code computes a solution to BL equation
%...
     close all
     clear all
%...
%... start a stopwatch timer
     tic
%...
%... set global variables
     global eps;
     global ne z0 zL z n method D1;
%...
%... spatial grid
     z0=0.0;
     zL=1.0;
     ne=1;
     n=201;
     dz=(zL-z0)/(n-1);
     z=[z0:dz:zL]';
%...
%... graphs of the flux function and derivative
     t = 0;
     f = flux(t,z);
     fx = dflux_dx(t,z);
     subplot(1,2,1)
     plot(z,f)
     ylabel('f(x)')
     xlabel('x')
     subplot(1,2,2)
     plot(z,fx)
     ylabel('f''(x)')
     xlabel('x')
%...
%... model parameter
%     eps = 0.0001;
     eps = 0.01;
%...
%... initial conditions
     x = 1 - 3*z;
     x(x<0) = 0;
%...
%... select slope limiter approximation
%...
%      method = 'Van_Leer'         %... Van Leer slope limiter
%      method = 'superbee'         %... Superbee slope limiter
%      method = 'smart'            %... Smart slope limiter
%      method = 'mc'               %... monotonized center limiter
%      method = 'minmod'           %... Minmod slope limiter
%      method = 'Koren'            %... Koren slope limiter
      method = 'KT'               %... Kurganov Tadmor centered scheme
%...
%... differentiation matrix (diffusive term - stagewise differentiation)
     D1 = five_point_centered_D1(z);
%...  
%... call to ODE solver 
%...
     t=[0:0.1:0.6];
     options = odeset('RelTol',1e-3,'AbsTol',1e-3,'Stats','on');
%...     
     [tout, xout] = ode15s('buckley_leverett_pde',t,x,options);
%...
%....plot results
     figure(2)
     plot(z,xout,'k');
     xlabel('z');
     ylabel('x(z,t)');
%     title('Buckley-Leverett equation')
     hold on
     axis([0 1 0 1.2]);
%...
     figure(3)
     surf(z,tout,xout);
     view(10,70);
     grid off;
     colormap winter
     axis([z0 zL 0 0.6 min(min(xout)) max(max(xout))]);
     xlabel('z');
     ylabel('t');
     zlabel('x(z,t)');
%...
%...
%... read the stopwatch timer
     tcpu=toc;

  
