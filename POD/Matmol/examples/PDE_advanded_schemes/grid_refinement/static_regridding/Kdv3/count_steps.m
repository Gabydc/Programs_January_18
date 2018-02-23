%... The MatMol Group (2016)
     function [value,isterminal,direction] = count_steps(t,x)
%...
     global maxsteps tprint tflag nsteps
%...
%... update step counter  
     nsteps = nsteps + 1;

     if nsteps < maxsteps || t < tprint
        value = 1;
     else
        value = 0;
        tflag = t - tprint;
     end    
     isterminal= 1;
     direction = 0;    
