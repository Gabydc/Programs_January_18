%... The MatMol Group (2016)
%...
%... Spruce budworm and forest... a problem with fast and slow dynamics
%...
%... The forest areas in Northern America are periodically defoliated 
%... by outbreaks of spruce budworm. The initiation of these outbreaks 
%... is controlled by the interactions between the slowly changing volume
%... of a growing forest susceptible to budworm, the more quickly changing
%... budworm densities and feeding responses of budworm's avian predators, 
%... and rapidly changing weather conditions.
%... These dyanmics were analysed in (Ludwig, Jones and Holling, 1978)
%...
%... the budworm's growth follows a logistic equation
%...
%... B  = rb*B*(1-B/Kb) - g
%...  t  
%...
%... where the carrying capacity
%...
%... Kb = k*S*E^2/(E^2+Te^2)
%...
%... is proportional to the amount of foliage available, i.e. proportional to
%... S, and depends on the physiological condition (energy) of the trees E 
%... (Kb declines sharply when E falls below a threshold Te)
%...
%... and the effect of predation is represented by
%...
%... g = beta*B^2/(alpha^2+B^2)
%...
%... this latter function has the following properties:
%...
%... a) the consumption of prey by individual predators (birds) 
%...    is limited by saturation to the level beta
%... b) there is a decrease in the effectiveness of predation at low
%...    prey density (here, like B^2), i.e., the birds have a variety 
%...    of alternative foods
%... c) alpha determines the scale of budworm densities at which saturation 
%...    takes place. This half-saturation density is proportional to the branch 
%...    surface area S, i.e. alpha = a*S 
%...
%... the total surface area of the branches in a stand takes the following form
%...
%... S  = rs*S*(1-(S/Ks)*(Ke/E))
%...  t  
%...
%... which allows S to approach its upper limit Ks. The additional factor Ke/E 
%... is inserted into the equation because S does not inevitably increase under 
%... stress conditions (death of branches or even whole trees)
%...
%... the energy reserve also satisfies an equation of the logistic type
%...
%... E  = re*E*(1-E/Ke) - P*B/S
%...  t  
%...
%... where the second term on the RHS describes the stress exerted on the trees 
%... by the budworm's consumption of foliage. In this expression B/S represents 
%... the number of budworms per branch. The proportionality factor P is given by
%... 
%... P = p*E^2/(Te^2+E^2)
%...
%... as the stress on the trees is related to the amont of foliage consumed 
%... (P declines sharply when E falls below a treshold Te)
%...
%... in these equations:
%...
%...  B       budworm density
%...
%...  S       total surface area of the branches in a stand
%...
%...  E       energy reserve
%...
%...  t       time
%...
%...  Kb      budworm's carrying capacity 
%...
%...  rb      budworm's specific growth rate
%...
%...  Ks      maximum branch density 
%...
%...  rs      specific branch growth rate
%...
%...  Ke      maximum energy level
%...
%...  re      specific E growth rate
%...
%...  alpha   half-saturation density
%...
%...  beta    maximum budworm predated
%...
%... the initial conditions (ICs) are taken as
%...
%... B(t=0) = 10    S(t=0) = 7000    E(t=0) = 1          
%...
%... Reference
%...
%... D. Ludwig, D.D. Jones and C.S. Holling,  
%... Qualitative Analysis of Insect Outbreak Systems: 
%... The Spruce Budworm and Forest 
%... Journal of Animal Ecology 47(1978), 315-332.
%...
%... the following code computes a solution to this problem
%...
     close all
     clear all
%...
%... start a stopwatch timer
     tic
%...
%... set global variables
     global rb k beta a rs Ks re Ke p Te
%...
%... model parameters
     rb = 1.52;
     k = 355;
     beta = 43200;
     a = 1.11;
     rs = 0.095; 
     Ks = 25440;
     re = 0.92;
     Ke = 1.0;
     p = 0.00195;
     Te = 0.03;
%...
%... initial conditions
     t0=0;
     tf=200;
     B = 10;
     S = 7000;
     E = 1;
     x=[B S E]';
%...  
%... call to ODE solver (comment/decomment one of the methods to select a solver)
%...
%      method = 'Euler'     
%...
%      method = 'rkf45'     
%...
%      method = 'ode45'     
%...
      method = 'ode15s'     
%...
     switch method
%...         
%... Euler
       case('Euler')
           Dt = 0.001;
           Dtplot = 0.5;
           [tout, yout] = euler_solver(@spruce_budworm_odes, t0, tf, x, Dt, Dtplot);
%... rkf45
       case('rkf45')
           hmin = 1e-3;
           nstepsmax = 10000;
           abstol = 1e-6;
           reltol = 1e-6;
           Dtplot = 0.5;
           [tout, yout, eout] = rkf45_solver(@spruce_budworm_odes,t0,tf,x,hmin,nstepsmax,abstol,reltol,Dtplot);
%... ode45
       case('ode45')
           options = odeset('RelTol',1e-6,'AbsTol',1e-6);
           t=[t0:0.5:tf];
           [tout, yout] = ode45(@spruce_budworm_odes,t,x,options);
%... ode15s
       case('ode15s')
           options = odeset('RelTol',1e-6,'AbsTol',1e-6);
           t=[t0:0.5:tf];
           [tout, yout] = ode15s(@spruce_budworm_odes,t,x,options);
%...
     end
%...
%... plot results
     figure(1)
     subplot(3,1,1)
     plot(tout,yout(:,1),'b');
%     xlabel('t [years]');
     ylabel('B(t)');
%     title('Budworm density');
     subplot(3,1,2)
     plot(tout,yout(:,2),'b');
%     xlabel('t [years]');
     ylabel('S(t)');
%     title('Branch density')
     subplot(3,1,3)
     plot(tout,yout(:,3),'b');
     xlabel('t [years]');
     ylabel('E(t)');
%     title('Energy');
     figure(2)
     plot3(yout(:,1),yout(:,2),yout(:,3),'b')
     xlabel('B(t)');
     ylabel('S(t)');
     zlabel('E(t)');
     grid
     title('3D phase plane plot');
%...
%... read the stopwatch timer
     tcpu=toc;


