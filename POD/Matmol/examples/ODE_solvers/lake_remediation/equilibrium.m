%... The MatMol Group (2016)
     function [N1star,N2star] = equilibrium(r0,Clim,alpha,Cdeath,K0,Klim,beta,H,Q,delta,delta0,eta,eta0)
%...
     Estar = eta*(Q/delta-Clim)/(eta0+eta*(delta0/delta));
%...
     Cstar = (Q-delta0*Estar)/delta
%...
     if Cstar < Clim
        r = r0;
     elseif Clim <= Cstar < Cdeath
        r = r0-alpha*(Cstar-Clim);
     else
        r = 0;
     end   
%...
     if Cstar < Clim
        K = K0;
     elseif Clim <= Cstar < Cdeath
        K = K0-beta*(Cstar-Clim);
     else
        K = Klim;
     end   
%...
     N1star  = (K*r)/(2*r0)*(1-sqrt(1-(4*r0*H)/(K*r^2)));
     N2star  = (K*r)/(2*r0)*(1+sqrt(1-(4*r0*H)/(K*r^2)));
