%... The MatMol Group (2016)
     function TL = temp(t)
%...
     global T1 T2 tlim Tad;

%...     
     if t<tlim
         TL = (T1+(T2-T1)*t/tlim)/Tad;
     else
         TL = T2/Tad;
     end
