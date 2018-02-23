%...  The MatMol Group (2016)
    function varargout = gni_fstrvl2(odefile,tspan,y0,options,varargin)
%GNI_FSTRVL2  Solve non-stiff differential equations.
%   [T,Q,P] = GNI_FSTRVL2('G',TSPAN,[Q0;P0]) with TSPAN = [T0 TFINAL] 
%   integrates the system of differential equations q'' = g(t,q) 
%   from time T0 to TFINAL with initial conditions q = Q0, q' = P0. 'G' is
%   a scalar. For a vector Y, G(~,Y) must return a column vector
%   corresponding to g(y). Each column in the solution array Y corresponds 
%   to a time returned in the column vector T.
%   
%   Use [T,Q] = GNI_FSTRVL2('G',TSPAN,[Q0;P0],OPTIONS) if you don't use P.
%
%   [T,Q,P] = GNI_FSTRVL2('G',TSPAN,[Q0;P0],OPTIONS) solves as above 
%   with default integration properties replaced by values in OPTIONS, 
%   an argument created with the GNI_SET function. Commonly used options 
%   are step size 'StepSize' (1e-2 by default) and method used 'Method'
%   ('817' by default).  
%
%   SOL = GNI_FSTRVL2(ODEFILE,TSPAN,Y0,...) returns a structure. The
%   structure SOL includes these fields:
%   sol.t          Steps chosen by the solver
%   sol.y          Each column sol.y(:,i) contains the solution at sol.t(i)
%   sol.solver     Solver name  
%   sol.stats      Statistics
%
%   Example   
%       options = gni_set('StepSize',0.01,'OutputFcn','odephas2','OutputSel',[1 2],'OutputSteps',10);
%       gni_fstrvl2('G_kepler_perturb',[],[],options);  

%     solves the system q'' = G_kepler_perturb(t,y), using step size defined
%     in the file 'G_kepler_perturb' and the default method '817' and 
%     plots the two first components of the solution in the phase space.
%
%   gni_fstrvl2 is an implementation of composition methods. It uses the
%   method of Störmer-Verlet as a basic integrator. It is for autonomous
%   second order Hamiltonian systems.
%   
%   REFERENCE:
%   E. Hairer, C. Lubich, G. Wanner, Geometric Numerical Integration. 
%   Structure-Preserving Algorithms for Ordinary Diferential Equations,
%   springer-verlag, berlin, second edition Edition, Vol. 31, Springer
%   Series in Computational Mathematics 31, 2006.

    solver_name = 'gni_fstrvl2';

% stats
    nfevals = 0; 
    
% Test that 'odefile' is a string
    if nargin == 0
        error('Not enough input arguments.  See gni_fstrvl2.');
    elseif ~ischar(odefile) && ~isa(odefile, 'inline')
        error('First argument must be a single-quoted string.  See gni_fstrvl2.');
    end

% Test that 'options' is a structure
    if nargin == 1
        tspan = []; y0 = []; options = [];
    elseif nargin == 2
        y0 = []; options = [];
    elseif nargin == 3
        options = [];
    elseif ~isempty(options) && ~isa(options,'struct')
        error('Correct syntax is gni_fstrvl2(''odefile'',tspan,y0,options).');
    end
    
% Output
    if nargout == 1
        sol = [];
        sol.solver = solver_name;
        sol.extdata.odefun = odefile;
        sol.extdata.options = options;
        sol.extdata.varargin = varargin;
    end
    
% Get default tspan and y0 from odefile if none are specified.
    if isempty(tspan) || isempty(y0)
        if (nargout(odefile) < 3) && (nargout(odefile) ~= -1)
            msg = sprintf('Use gni_fstrvl2(''%s'',tspan,y0,...) instead.',odefile);
            error(['No default parameters in ' upper(odefile) '.  ' msg]);
        end
        [def_tspan,def_y0,def_options] = feval(odefile,[],[],'init',varargin{:});
        nfevals = nfevals + 1;
        if isempty(tspan)
            tspan = def_tspan;
        end
        if isempty(y0)
            y0 = def_y0;
        end
        if isempty(options)
            options = def_options;
        else
            options = gni_set(def_options,options);
        end
    end

% Test that tspan is internally consistent.
    ntspan = length(tspan);
    if ntspan > 2
        error('The timespan must consist of one or two numbers.');
    end
    if ntspan == 1
        tspan = [0 tspan];
        t0 = 0;
    else
        t0 = tspan(1);
    end
    tf = tspan(2);
    if t0 >= tf
        error('The final time must greater than the starting time.');
    end

% Test that the length of y0 is pair
    if size(y0,2) ~= 1
        y0 = y0';
    end
    ny = length(y0);
    if (mod(ny,2) ~= 0)
        error('The initial value must have 2*N components, where N is the dimension of the system.');
    end
    npq = ny / 2;
    
% Get parameters
    h = gni_get(options,'StepSize');
    hold = h;
    if (isempty(h))
        nsteps = gni_get(options,'NumSteps');
        if (isempty(nsteps))
            h = 1e-2;
            fprintf('Warning: No initial step size provided, using h = %e instead.\n',h);
        else
            h = (tf - t0)/nsteps;
        end
    end
    if ((~isa(h,'double')) || (size(h,1) ~= 1) || (size(h,2) ~= 1)) 
        error('The option ''Precision'' must contain a single number');
    end
    nt = round((tf-t0)/h) + 1;
    h = (tf-t0) / (nt-1);  
    if hold~=h
        fprintf('Integration period isn''t a multiple integer of the StepSize, \n')
        fprintf('take h=%g instead \n',h)
    end
    
% Test that odefile is consistent.
    Qn = y0(1:npq);
    Pn = y0(npq+1:ny);
    if isempty(varargin)
        args = {};
    else
        args = [{[]} varargin];
    end
    F0 = feval(odefile,[],Qn,args{:});
    nfevals = nfevals + 1;
    [mf,nf] = size(F0);
    if nf > 1
        error([upper(odefile) ' must return a column vector.'])
    elseif mf ~= npq
        msg = sprintf('an initial condition vector of length 2*%d.',mf);
        error(['Solving ' upper(odefile) ' requires ' msg]);
    end

% Outputstep.
    outstep = gni_get(options,'OutputSteps',1);
    if ((~isa(outstep,'double')) || size(outstep,1) ~= 1 || size(outstep,2) ~= 1)
        error('The option ''OutputSteps'' must contain a single integer');
    end
    if (outstep <= 0)
        outstep = 1;
    end
    outstep = round(outstep);    
    
% If we are generating output.
    if nargout > 0
        nout = 1;
        
        Tout = t0:outstep*h:tf;
        Qout = zeros(npq,length(Tout));
        Pout = zeros(npq,length(Tout));
        Qout(:,1) = Qn; 
        Pout(:,1) = Pn; 
    end
    
% Ouput fuction.
    if nargout > 0
        haveoutput = true;
        outfun = gni_get(options,'OutputFcn');
    else
        haveoutput = false;
        outfun = gni_get(options,'OutputFcn','odeplot');
    end
    if isempty(outfun)
        haveoutfun = false;
    else
        haveoutfun = true;
        outputs = gni_get(options,'OutputSel',1:ny);
        feval(outfun,[t0 tf],y0(outputs),'init');
    end       
    
% Initialize the method    
    method = gni_get(options,'Method','817');
    gamma = coeff_comp(method);
    ng = max(size(gamma));
    hg = h * gamma;
    h2 = hg / 2;   
    Ep = 0;
    Eq = 0;
    
%% THE MAIN LOOP
    outpoint = 0;
    for i = 1:nt-1
        addpoint = false;
        outpoint = outpoint + 1;
        if (outpoint >= outstep)
            outpoint = 0;
            addpoint = true;
        end        
        
        for j = 1:ng  
            Eq = Eq + h2(j) * Pn;
            q12 = Qn + Eq;
            Eq = Eq + (Qn - q12);

            Pold = Pn;
            Ep = Ep + hg(j) * feval(odefile,[],q12);
            Pn = Pold + Ep;
            Ep = Ep + (Pold - Pn);

            Eq = Eq + h2(j) * Pn;          
            Qn = q12 + Eq;
            Eq = Eq + (q12 - Qn);
        end
        
    % If we are generating output.
        if addpoint
            if haveoutput
                nout = nout + 1;
                Qout(:,nout) = Qn;
                Pout(:,nout) = Pn;
            end
            if haveoutfun
                t = t0 + i*h;
                Yn = [Qn;Pn];
                if feval(outfun,t,Yn(outputs),'') == 1
                    return;
                end
            end
        end
    end  
    
%% Problem solved

if haveoutfun
  feval(outfun,[],[],'done');
end

% Finalize the output
    solver_output = {};
    if nargout == 1
        sol.t = Tout;
        sol.q = Qout;
        sol.p = Pout;  
        sol.stats.nfevals = nfevals + (nt-1) * ng;
        sol.stats.nsteps = nt;
        solver_output{1} = sol;
    elseif nargout==2
        solver_output{1} = Tout;
        solver_output{2} = Qout;
    elseif nargout == 3
        solver_output{1} = Tout;
        solver_output{2} = Qout;
        solver_output{3} = Pout;       
    end

    varargout = solver_output;
    