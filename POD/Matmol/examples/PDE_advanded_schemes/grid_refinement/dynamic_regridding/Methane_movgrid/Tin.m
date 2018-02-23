%... The MatMol Group (2016)
     function TL = Tin(t)
%...
     tau = 10^(-2);
%...     
     if t < 10
         TL = 300;
     else
         TL = 300+200*(1-exp(-(t-10)/tau));
     end
