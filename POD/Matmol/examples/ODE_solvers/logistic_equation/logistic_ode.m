%... The MatMol Group (2016)
     function Nt = logistic_ode(t,N)
%...
%... Set global variables
     global a b K N0
%...
     Nt  = (a-b*N)*N;
