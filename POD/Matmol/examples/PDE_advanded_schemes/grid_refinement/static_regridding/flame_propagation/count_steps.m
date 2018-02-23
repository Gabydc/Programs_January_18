%... The MatMol Group (2016)
     function [value,isterminal,direction] = count_steps(t,x)
%...
     global nsteps maxsteps
%...
%... update step counter
     nsteps = nsteps + 1;
%...
%... check if number of steps exceeds maxsteps
     value = maxsteps-nsteps;
     isterminal = 1;    
     direction = 0;
     if value == 0,   
         nsteps = 0;
     end    