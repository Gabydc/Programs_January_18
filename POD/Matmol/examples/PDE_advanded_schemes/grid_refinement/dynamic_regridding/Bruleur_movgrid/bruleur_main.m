%... The MatMol Group (2016)
%...
%... Combustion process of Aluminium dust
%...
%...          
%... c  = -S(c,T)
%...  t     
%...
%... T  = +S(c,T) + psi*T
%...  t                  zz
%...
%... with S(c,T)=2*c^2*exp[gam*(1-1/T)]
%...
%... where
%...
%...  c  : Al concentration : dependent variable
%...  T  : temperature : dependent variable
%...
%...  t       time
%...
%...  z       space : 0 < z < 2
%...
%... I.C. : c(z,0)=1 ;  T(z,0)=T1/Tad
%...                                                                  t     1
%... B.C. : T(0,t)=temp(t)  ;  T (L,t)=0  with  temp(t)= [T1+(T2-T1)*----]*---  if  t < tlim
%...                            z                                    tlim  Tad
%...                                                      T2
%...                                                   = ---   if   t > tlim
%...                                                     Tad
%...
%... numerical values  :
%...
%...    gam = 11.35
%...
%...    psi = 1.89*10e-4
%...
%...    T1 = 725
%...
%...    T2 = 1275
%...
%...    Tad = 1854.89
%...
%...    tlim=10e-3
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
     global psi T1 T2 tlim Tad gam
%...
%... temporal conditions
%...
     t=[0:10:200];
%...
%... moving grid parameters
     alpha=1;
     kappa=1; 
     mu=kappa*(kappa+1);
     tau=1e-2;
%...
%... PDE parameters
     ne=2;
     gam=1.75*90/(8.314*1.6694);
     psi=10^(-8)/(4.5*exp(-gam));
     T1=725;
     T2=1275;
     Tad=16694/9;
     tlim=10^(-2);
%...
%... initial spatial grid z : if n is the number of moving grid points,
%... dim(z) = n+2
%...
     zL=0.0;
     zR=2.0;
     n=99;
     dz=(zR-zL)/(n+1);
     z=[zL:dz:zR]';
%...
%... initial dependent variables u(n+2,ne)
     u(1:n+2,1)=1;
     u(1:n+2,2)=T1/Tad;
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
%... 2d spatial derivative by 3 pts centered finite difference
     eq(:,:,1) = [0 0 1;0 0 1];
     eq(:,:,2) = [0 0 1;-1 1 1];
%...
%... monitoring choice : 
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
         ordo1(i,1)= ordo1(i,2)-(absc(i,2)-absc(i,1))*(ordo1(i,3)-ordo1(i,2))/(absc(i,3)-absc(i,2));
         ordo1(i,n+2)=ordo1(i,n+1);
         ordo2(i,2:n+1)=Tad*yout(i,ne+1-ne+1:ne+1:n*(ne+1)-ne+1);
         ordo2(i,1)=Tad*temp(tout(i));
         ordo2(i,n+2)=ordo2(i,n+1);
     end

     axes(h1);
     hold on
     for i=1:length(tout),
         plot(absc(i,:),ordo1(i,:),'.-k')
     end
     xlabel('z');
     ylabel('c(z,t)');
%     title('combustion de poudre de Al')
     axis([0 2 -.1 1.1])
 
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
     ylabel('T(z,t)');
%     title('combustion de poudre de Al')

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
