%... The MatMol Group (2016)
     function xt = spring_mass_odes(t,x)
%...
%... Set global variables
     global m k
%...
%... Transfer dependent variables  
     z = x(1);
     v = x(2);
%...      
%... Temporal derivatives
%...
     zt  = v;
     vt  = -(k/m)*z;
%...
%... Transfer temporal derivatives
     xt = [zt vt]';
