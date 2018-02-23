%... The MatMol Group (2016)
%... Equillibrium point for the spruce budworm problem
     clear all

     clear all

%... start a stopwatch timer
     tic

%... set global variables
     global rb k beta a

%... model parameters
     rb    = 1.52;
     k     = 355;
     beta  = 43200;
     alpha = 10000;

%... plot Bt = f(B) for various values of S
     B = [0:1000:200000];
     for S = 200:200:10000
    %... Computation of dB/dt
         Kb = k*S;
         Bt = rb*B.*(1-B/Kb) - (beta*B.^2)./(alpha^2+B.^2);
         figure (1)
         plot(B,Bt);
         hold on
         plot(B,zeros(size(B)),'-k')
    
         axis([0 2e5 -1e5 1e5]);
         xlabel('B','FontName','Helvetica','FontSize',12)
         ylabel('dB/dt','FontName','Helvetica','FontSize',12)
    
    %... Compute the roots of the equation
         a = rb;
         b = Kb;
         c = beta;
         d = alpha;
         roots = solve('a*x*(1-x/b) - (c*x^2)/(d^2+x^2) = 0');
         r(1) = eval(roots(1));
         r(2) = eval(roots(2));
         r(3) = eval(roots(3));
         r(4) = eval(roots(4));
         plot(r,zeros(size(r)),'ok')
     end    
