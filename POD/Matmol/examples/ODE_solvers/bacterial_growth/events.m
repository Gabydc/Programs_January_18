%... The MatMol Group (2016)
    function [value,isterminal,direction] = events(t,x)     
%...
%... Transfer dependent variables
     X = x(1);
     S = x(2);
%...
%... check if substrate concentration becomes negative
     value = S;         % monitor S and see if it vanishes
     isterminal = 1;    % stop integration when S vanishes
     direction = -1;    % S decreases from initial positive values