%... The MatMol Group (2016)
%...
%...          
%... c  = Deff * c   - v0 * c  - r(c,T)
%...  t           zz        z
%...
%...       lamb          eps*v0*rhog*cpg         2*kw                -DH
%... T  = ------ * T   - --------------- * T  + --------- * (Tw-T) + ------ * r(c,T)
%...  t   rho*cp    zz        rho*cp        z   Rt*rho*cp            rho*cp
%...
%...                    c*exp[-Eact/(R*T)]
%... with r(c,T) = Kr * ------------------
%...                       1 + Kc*c
%...
%... where
%...
%...  c      : CO2 concentration : dependent variable
%...
%...  T      : temperature : dependent variable
%...
%...  t      : time
%...
%...  z      : space : 0 < z < .2
%...
%...  Deff   : effective diffusion coefficient : 5*10e-4
%...
%...  rho*cp : capacit� calorifique du r�acteur : 364 (not� rhocp ci-dessous)
%...
%...  rhog   : masse volumique du gaz : 0.0775
%...
%...  cpg    : chaleur sp�cifique du gaz : 2.29
%...
%...  Tw     : T� de la paroi du r�acteur : 300
%... 
%...  Eact   : �nergie d'activation : 25211
%...
%...  -DH    : chaleur de r�action : 6006 (not� DHM ci-dessous)
%...
%...  lamb   : conductibilit� : 3.5*10-4
%...
%...  v0     : vitesse du gaz : 1.5
%...
%...  eps    : porosit� du substrat : 0.6
%...
%...  kw     : coefficient d'�change de chaleur : 5*10e-4
%...
%...  Rt     : rayon du reacteur (tubulaire) : 0.01
%...
%...  Kr     : gain du terme de r�action : 0.971*10e13
%...
%...  Kc     : autre coeff. du terme de r�action : 12.7
%...
%...  R      : cte des gaz parfaits : 1.98
%...
%...  I.C.   : c(z,0)=0 ;  T(z,0)=300
%...
%...                              v                       eps*v*rhog*cpg
%...  B.C.   : at z = 0  : c  = ---- * (c-cin)   ;   T  = -------------- * (T-Tin)
%...                        z   Deff                  z        lamb
%...
%...           at z = .2 : c  = T  = 0
%...                        z    z
%...
%...           where cin = 2.5*u(t-10) and  Tin=300+200*u(t-10)
%...
%...
%... Implementation of MOVGRID - Finite differences version
%...
%...
     close all
     clear all
%...
%... start a stopwatch timer
%...
     tic
%... 
%... set global variables
%...
     global zL zR n ne X eq choice
     global alpha kappa mu tau
     global u z
     global Deff rhocp rhog cpg Tw Eact DHM lamb v0 eps kw Rt Kr Kc R;
%...
%... temporal conditions
     t=[0:50:1000];
     dt=50;
%...
%... moving grid parameters
     alpha=50;
     kappa=2; 
     mu=kappa*(kappa+1);
     tau=1e-2;
%...
%... PDE parameters
     ne=2;
     Deff = 5*(10^(-4));
     rhocp = 364;
     rhog = 0.0775;
     cpg = 2.29;
     Tw = 300;
     Eact = 25211;
     DHM = 6.006;
     lamb = 3.5*(10^(-4));
     v0 = 1.5;
     eps = 0.6;
     kw = 5*(10^(-4));
     Rt = 0.01;
     Kr = 0.971*(10^13);
     Kc = 12.7;
     R = 1.98;
%...
%... initial spatial grid z : if n is the number of moving grid points,
%... dim(z) = n+2
     zL=0.0;
     zR=.2;
     n=31;
     dz=(zR-zL)/(n+1);
     z=[zL:dz:zR]';
%...
%... initial dependent variables u(n+2,ne)
     u(1:n+2,1)=0;
     u(1:n+2,2)=300;
%...
%... global variables for JPat : function JPat implements the sparsity
%... pattern of the jacobian in X ; informations depending on the
%... choosen spatial derivative matrices are stocked in the ne eq(ne,ne,3) matrices.
%... 
%...    example : hypothetic problem : 
%...
%...        u1  = -f1  + eps1*u1   + k1*u2      where  f1 depends on u3
%...          t      z          zz
%...
%...        u2  = -f2  + eps2*u2                where  f2 depends on u2 and u3
%...          t      z          zz
%...
%...        u3  = -f3  + eps3*u3                where  f3 depends on u1, u2 and u3
%...          t      z          zz
%...
%...    where the 1st (2nd) derivative is computed by a 3 (5) pts centered finite difference
%... 
%...    It results :
%...
%...        eq(:,:,1) = [-2 2 1; 0 0 1;-1 1 1];
%...        eq(:,:,2) = [ 0 0 0;-2 2 1;-1 1 1];
%...        eq(:,:,3) = [-1 1 1;-1 1 1;-2 2 1];

     X=zeros(n*(ne+1),n*(ne+1));
%... 
%... see right_vect file :
%... convection term : 5 pts biased upwind finite difference
%... diffusion term : 5 pts centered finite difference
     eq(:,:,1) = [-3 2 1;0 0 1];
     eq(:,:,2) = [0 0 1;-3 2 1];
%...
%... monitoring choice
%...
%...                                                     1    ne (u(i+1,j)-u(i,j))^2
%...    choose between : 'der1'  : m(i) = sqrt[ alpha + ---  sum ------------------- ]
%...                                                     ne  j=1   (z(i+1)-z(i))^2
%...
%...        
%...
%...                                                     1    ne     uzz(i+1,j)+uzz(i,j)
%...                     'der2'  : m(i) = sqrt[ alpha + ---  sum abs(-------------------)]
%...                                                     ne  j=1           2  
%...
%...
%...
%...                                                     1    ne (fl[u(i+1,j)]-fl[u(i,j)])^2
%...                     'derfl' : m(i) = sqrt[ alpha + ---  sum --------------------------- ]
%...                                                     ne  j=1      (z(i+1)-z(i))^2
%...

     choice = 'der1';
%...
%... global initial condition
%...
     for j=1:ne,
         x(j:ne+1:n*(ne+1))=u(2:n+1,j);
     end
     x(ne+1:ne+1:n*(ne+1))=z(2:n+1);
%...
%... %%%%%%%%%%%%%%%%%%%%%%
%... % call to ODE solver %
%... %%%%%%%%%%%%%%%%%%%%%%
%...
     options = odeset('RelTol',1e-6,'AbsTol',1e-6,'Mass',@mass,...
         'MStateDependence','strong','JPattern',JPat,'MvPattern',MvPat,'MassSingular','no','stats','on');
%...     
     [tout, yout] = ode15s(@right_vect,t,x,options);
%...
%... %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%... % visualisation of the solution %
%... %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%...
     h1=axes('position',[0.1 0.22 0.8 0.75]);
     h2=axes('position',[0.1 0.1 0.8 0.10]);

     nt=length(tout);
     absc=zeros(nt,n+2);
     ordo1=zeros(nt,n+2);
     ordo2=zeros(nt,n+2);

     for i=1:nt
         absc(i,2:n+1)=yout(i,ne+1:ne+1:n*(ne+1));
         absc(i,1)=zL;
         absc(i,n+2)=zR;
         ordo1(i,2:n+1)=yout(i,ne+1-ne:ne+1:n*(ne+1)-ne);
         ordo1(i,1)=(ordo1(i,2)+v0*(absc(i,2)-absc(i,1))*cin((i-1)*dt)/Deff)/(1+v0*(absc(i,2)-absc(i,1))/Deff);
         ordo1(i,n+2)=ordo1(i,n+1)+(absc(i,n+2)-absc(i,n+1))*(ordo1(i,n+1)-ordo1(i,n))/(absc(i,n+1)-absc(i,n));
         ordo2(i,2:n+1)=yout(i,ne+1-ne+1:ne+1:n*(ne+1)-ne+1);
         ordo2(i,1)=(ordo2(i,2)+eps*v0*rhog*cpg*(absc(i,2)-absc(i,1))*Tin((i-1)*dt)/lamb)/...
                    (1+eps*v0*rhog*cpg*(absc(i,2)-absc(i,1))/lamb);
         ordo2(i,n+2)=ordo2(i,n+1)+(absc(i,n+2)-absc(i,n+1))*(ordo2(i,n+1)-ordo2(i,n))/(absc(i,n+1)-absc(i,n));
     end

     axes(h1);
     hold on
     for i=1:length(tout),
         plot(absc(i,:),ordo1(i,:),'.-k')
     end
     xlabel('z');
     ylabel('c');

     axes(h2);
     hold on
     for i=1:length(tout),
         plot(absc(i,:),((i-1)*1.)*ones(1,n+2),'.k')
         hold on
     end
     xlabel('z');
     ylabel('t');
     
     figure
     h1=axes('position',[0.1 0.22 0.8 0.75]);
     h2=axes('position',[0.1 0.1 0.8 0.10]);
     axes(h1);
     hold on
     for i=1:length(tout),
         plot(absc(i,:),ordo2(i,:),'.-k')
     end
     xlabel('z');
     ylabel('T');

     axes(h2);
     hold on
     for i=1:length(tout),
         plot(absc(i,:),((i-1)*1.)*ones(1,n+2),'.k')
         hold on
     end
     xlabel('z');
     ylabel('t');
%...     
%... read the stopwatch timer
%...
     tcpu=toc    
%... 
