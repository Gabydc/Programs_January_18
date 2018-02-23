%...  The MatMol Group (2016)
     function Stats = AddIntStats(Stats1,Stats2,solver)

%... Returns a structure Stats which contains the sum of the ODE integration 
%... statistics for two structures Stats1 and Stats2 containing the individual
%... integration statistics.

    switch lower(solver)
        case {'ode113','ode23','ode45'}
            Stats.nsteps   = Stats1.nsteps   + Stats2.nsteps;
            Stats.nfailed  = Stats1.nfailed  + Stats2.nfailed; 
            Stats.nfevals  = Stats1.nfevals  + Stats2.nfevals;        
        case {'ode15i','ode15s','ode23s','ode23t','ode23tb'}
            Stats.nsteps   = Stats1.nsteps   + Stats2.nsteps;
            Stats.nfailed  = Stats1.nfailed  + Stats2.nfailed; 
            Stats.nfevals  = Stats1.nfevals  + Stats2.nfevals; 
            Stats.npds     = Stats1.npds     + Stats2.npds;
            Stats.ndecomps = Stats1.ndecomps + Stats2.ndecomps;
            Stats.nsolves  = Stats1.nsolves  + Stats2.nsolves;
        otherwise
            disp('Wrong ODE solver.')
    end
