     function [value,isterminal,direction] = events_StatRegrid(t,x)
%...
%... The MatMol Group (2009)
%...
     global nsteps maxsteps
%...
%... update step counter  
     nsteps = nsteps + 1;
%...
%... check if maximum number of steps (maxsteps) is exceeded
     if nsteps <= maxsteps
        value = 1;
     else
        value = 0;
     end    
     isterminal = 1;    
     direction = 0;    
