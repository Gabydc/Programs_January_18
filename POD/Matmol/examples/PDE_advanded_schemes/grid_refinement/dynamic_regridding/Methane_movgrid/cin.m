%... The MatMol Group (2016)
     function CL = cin(t)
%...
     tau = 10^(-2);
%...     
     if t < 10
         CL = 0;
     else
         CL = 2.5*(1-exp(-(t-10)/tau));
     end
