%...  The MatMol Group (2016)
    function varargout = gni_vstrvl2(odefile,tspan,y0,options,varargin)

%GNI_VSTRVL2  Solve non-stiff differential equations.
%   [T,Q,P] = GNI_VSTRVL2('G',TSPAN,[Q0;P0]) with TSPAN = [T0 TFINAL] 
%   integrates the system of differential equations q'' = g(t,q) 
%   from time T0 to TFINAL with initial conditions q = Q0, q' = P0. 'G' is
%   a scalar. For a vector Y, G(~,Y) must return a column vector
%   corresponding to g(y). Each column in the solution array Y corresponds 
%   to a time returned in the column vector T.
%   
%   Use [T,Q] = GNI_VSTRVL2('G',TSPAN,[Q0;P0],OPTIONS) if you don't use P.
%
%   [T,Q,P] = GNI_VSTRVL2('G',TSPAN,[Q0;P0],OPTIONS) solves as above 
%   with default integration properties replaced by values in OPTIONS, 
%   an argument created with the GNI_SET function. Commonly used options 
%   are step size 'Precision' (1e-1 by default) and method used 'Method'
%   ('817' by default). You can also change de time rescaling function with
%   'RescalFun'.
%
%   SOL = GNI_VSTRVL2(ODEFILE,TSPAN,Y0,...) returns a structure. The
%   structure SOL includes these fields:
%   sol.t          Steps chosen by the solver
%   sol.y          Each column sol.y(:,i) contains the solution at sol.t(i)
%   sol.solver     Solver name  
%   sol.stats      Statistics
%
%   Examples:
%   1)  gni_vstrvl2('G_kepler_perturb')
%
%   2)  options = gni_set('RescalFun','resc_lenard');
%       gni_vstrvl2('G_lennard',[],[],options);  
%
%     solves the system q'' = G_lennard(t,y), using step size function 
%     defined in the file 'resc_lenard' and the default method '817' and 
%     plots Q and P as a function of time.
%
%   gni_vstrvl2 is an implementation of composition methods in variable 
%   step size. It uses the method of Störmer-Verlet as a basic integrator. 
%   It is for autonomous second order Hamiltonian systems.
%   
%   REFERENCES:
%   E. Hairer, C. Lubich, G. Wanner, Geometric Numerical Integration. 
%   Structure-Preserving Algorithms for Ordinary Di?erential Equations,
%   springer-verlag, berlin, second edition Edition, Vol. 31, Springer
%   Series in Computational Mathematics 31, 2006.

solver_name = 'gni_vstrvl2';

% stats
    nfevals = 0;

% Test that 'odefile' is a string.
    if nargin == 0
        error('Not enough input arguments.  See gni_vstrvl2.');
    elseif ~ischar(odefile) && ~isa(odefile, 'inline')
        error('First argument must be a single-quoted string.  See gni_vstrvl2.');
    end

% Test that 'options' is a structure.
    if nargin == 1
        tspan = []; y0 = []; options = [];
    elseif nargin == 2
        y0 = []; options = [];
    elseif nargin == 3
        options = [];
    elseif ~isempty(options) && ~isa(options,'struct')
        error('Correct syntax is gni_vstrvl2(''odefile'',tspan,y0,options).');
    end
    
% Output.
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
            msg = sprintf('Use gni_vstrvl2(''%s'',tspan,y0,...) instead.',odefile);
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
    end
    if tspan(1) >= tspan(2)
        error('The final time must greater than the starting time.');
    end    
    t = tspan(1);
    tfin = tspan(2);

% Test that the length of y0 is pair.
    if min(size(y0)) ~=1
        error('The initial value must be a vector column.');
    end
    if size(y0,2) ~= 1 % verifier qu'une matrice ne peut pas passer...
        y0 = y0';
    end
    ny = length(y0);
    if (mod(ny,2) ~= 0)
        error('The initial value must have 2*N components, where N is the dimension of the system.');
    end
    
% get parameters.
    prec = gni_get(options,'Precision');
    if (isempty(prec))
            prec = 1e-1; 
            fprintf('Warning: No precision parameter provided, using prec = %e instead.\n',prec);
    end
    if ((~isa(prec,'double')) || (size(prec,1) ~= 1) || (size(prec,2) ~= 1)) 
        error('The option ''Precision'' must contain a single number');
    end
    
% Test that odefile is consistent.  
    npq = ny / 2;
    Qn = y0(1:npq);
    Pn = y0(npq+1:ny);
    if isempty(varargin)
        args = {};
    else
        args = [{[]} varargin];
    end
    F1 = feval(odefile,[],Qn,args{:});
    nfevals = nfevals + 1;
    [mf,nf] = size(F1);
    if nf > 1
        error([upper(odefile) ' must return a column vector.'])
    elseif mf ~= npq
        msg = sprintf('an initial condition vector of length 2*%d.',mf);
        error(['Solving ' upper(odefile) ' requires ' msg]);
    end
    
% outputstep
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
        
        Qout = zeros(npq,1);
        Pout = zeros(npq,1);
        Qout(:,1) = Qn;
        Pout(:,1) = Pn;
    
        Tout = t;
    end
    
% Ouput fuction 
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
        feval(outfun,[t tfin],y0(outputs),'init');
    end    

% Initialize the method.    
    method = gni_get(options,'Method','817');
    gamma = coeff_comp(method);
    ng = size(gamma);
    ng = max(ng);
    EQ = 0;
    EP = 0;      
    
% Initialise du rescaling function.
    RescalFun = gni_get(options,'RescalFun');
    if (isempty(RescalFun)) 
        haveRescalFun = false;
    else
        haveRescalFun = true;
    end   
    
% Compute s12 for first step.
    if haveRescalFun 
        s12 = feval(RescalFun,t,Qn,Pn,args{:});
    else
        s12 = norm([Pn;F1]); 
    end
    
%% THE MAIN LOOP.
    outpoint = 0;
    n = 1;
    done = true;
    while done
        told = t;
        
        addpoint = false;
        outpoint = outpoint + 1;
        if (outpoint >= outstep)
            outpoint = 0;
            addpoint = true;
        end
        
    % Compute step size.
        h = prec / s12;
        t = told + h;
        
    % Check final time.
        if t >= tfin            
            if t == tfin
                break;
            else            
                h = tfin - told;     % correction if overflow
                t = tfin;
                done = false;
            end            
        end
        
    % Start of actual step ynew = yold + h*...
        hg = h * gamma;
        hg2 = hg / 2;
        for i = 1:ng  
            EQ = EQ + hg2(i) * Pn;
            Q12 = Qn + EQ;
            EQ = EQ + (Qn - Q12);

            Pold = Pn;
            EP = EP + hg(i) * feval(odefile,[],Q12,args{:});
            Pn = Pold + EP;
            EP = EP + (Pold - Pn);

            EQ = EQ + hg2(i) * Pn;          
            Qn = Q12 + EQ;
            EQ = EQ + (Q12 - Qn);
        end  
        nfevals = nfevals + ng;
     
     % If we generate output
        if addpoint
            if haveoutput
                nout = nout + 1;

                Qout(:,nout) = Qn;
                Pout(:,nout) = Pn;   

                Tout(nout) = t;
            end
            if haveoutfun
                Yn = [Qn;Pn];
                if feval(outfun,t,Yn(outputs),'') == 1
                    return;
                end
            end
        end
        
    % Compute s12 for next step.
        if haveRescalFun 
            sig_inv = feval(RescalFun,t,Qn,Pn,args{:}); % must return the inverse of sig
        else
            F1 = feval(odefile,[],Qn,args{:});
            sig_inv = norm([Pn;F1]);
            nfevals = nfevals + 1;
        end       
        s12 = 2 * sig_inv - s12;  
        
        n = n + 1;
    end 
    
%% Problem solved.

if haveoutfun
    feval(outfun,[],[],'done');
end

% Finalize the output.
    solver_output = {};
    if nargout == 1
        sol.t = Tout;
        sol.q = Qout;
        sol.p = Pout;        
        sol.stats.nfevals = nfevals;
        sol.stats.nsteps = n;
        solver_output{1} = sol;
    elseif nargout == 2
        solver_output{1} = Tout;
        solver_output{2} = Qout;
    elseif nargout == 3
        solver_output{1} = Tout;
        solver_output{2} = Qout;
        solver_output{3} = Pout;
    end
    
    varargout = solver_output;
    