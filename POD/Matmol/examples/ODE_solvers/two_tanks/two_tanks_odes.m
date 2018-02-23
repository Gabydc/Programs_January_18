%... The MatMol Group (2016)
     function xt = Two_tanks_odes(t,x)     
%...
%... set global variables
     global A B eps K
%...      
%... temporal derivatives
%...
     u   = -K*x(1:2);
     xt  = A*x+B*u;
