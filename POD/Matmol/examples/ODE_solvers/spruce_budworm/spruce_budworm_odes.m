%... The MatMol Group (2016)
     function xt = spruce_budworm_odes(t,x)
%...
%... Set global variables
     global rb k beta a rs Ks re Ke p Te
%...
%... Transfer dependent variables  
     B=x(1);
     S=x(2);
     E=x(3);
%...      
%... Temporal derivatives
%...
     Kb = k*S*E^2/(E^2+Te^2);
     alpha = a*S;
     g = beta*B^2/(alpha^2+B^2);
     P = p*E^2/(Te^2+E^2);
%...
     Bt  = rb*B*(1-B/Kb) - g;
     St  = rs*S*(1-(S/Ks)*(Ke/E));
     Et  = re*E*(1-E/Ke) - P*B/S;
%...
%... Transfer temporal derivatives
     xt=[Bt St Et]';
