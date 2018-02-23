%... The MatMol Group (2016)
%... Program to compute the characteristic curves of the advection
%... equation
     clear all
     clc

%... Advection equation with constant velocity v=2
%... The equation of the characteristic equation is:
%... z = v*t + xi,   or   t = 1/v*(z-xi) 
%... where xi is how much the characteristic curve is propagated
     v    = 2;
     z_v2 = 0:0.1:1;
     xi   = 0;
     for ii=1:10
         t_v2(ii,:) = 1/v*(z_v2-xi);
         xi      = xi+0.1;
     end

     figure
     plot(z_v2,t_v2,'b')
     axis([0 1 0 0.6])

%... Advection equation with velocity v=2*z
%... The equation of the characteristic equation is:
%... z = xi*exp(beta*t),   or   t = 1/beta*ln(z/xi)
%... where xi is how much the characteristic curve is propagated
     z_vz = 0:0.02:1;
     xi   = 0.1;
     beta = 2;
     for ii=1:10
         t_vz(ii,:) = 1/beta*log(z_vz/xi);
         xi         = xi+0.1;
     end

     figure
     plot(z_vz,t_vz,'b')
     axis([0 1 0 1.2])
