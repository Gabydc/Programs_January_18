%...  The MatMol Group (2016)
    function varargout = gni_fgauss1(odefile,tspan,y0,options,varargin)
%GNI_FGAUSS1  Solve non-stiff differential equations.
%   [T,Y] = GNI_FGAUSS1('F',TSPAN,Y0) with TSPAN = [T0 TFINAL] integrates 
%   the system of differential equations y' = f(t,y) from time T0 to TFINAL
%   with initial conditions Y0. 'F' is a string. For a scalar t and a 
%   vector y, f(t,y) must return a column vector. Each column in the
%   solution array Y corresponds to a time returned in the column vector T.
%   
%   [T,Y] = GNI_FGAUSS1('F',TSPAN,Y0,OPTIONS) solves as above 
%   with default integration properties replaced by values in OPTIONS, 
%   an argument created with the GNI_SET function. Commonly used options 
%   are step size 'StepSize' (1e-2 by default) and method used 'Method'
%   ('G12' by default).  
%
%   SOL = GNI_FGAUSS1('F',TSPAN,Y0,...) returns a structure. The
%   structure SOL includes these fields:
%   sol.t          Steps chosen by the solver
%   sol.y          Each column sol.y(:,i) contains the solution at sol.t(i)
%   sol.solver     Solver name  
%   sol.stats      Statistics
%
%   Example    
%         [t,y]=gni_fgauss1('F_lotka',[0 20],[0.7 0.7]);   
%         plot(t,y(1,:));
%     solves the system y' = F_lotka(t,y), using the default step size 
%     1e-2 and the method 'G12' and plots the first component of the 
%     solution. 
%
%   gni_fgauss1 is an implementation of the Gauss methods. It uses fixed
%   point iteration to solve the nonlinear system at every step. It is
%   for non-stiff Hamiltonian and/or reversible systems.
%   
%   REFERENCE:
%   E. Hairer, C. Lubich, G. Wanner, Geometric Numerical Integration. 
%   Structure-Preserving Algorithms for Ordinary Di?erential Equations,
%   springer-verlag, berlin, second edition Edition, Vol. 31, Springer
%   Series in Computational Mathematics 31, 2006.

solver_name = 'gni_fgauss1';

% Stats.
    nfevals = 0; 
    niters = 0;
    
% Test that 'odefile' is a string.
    if nargin == 0
        error('Not enough input arguments.  See gni_fgauss1.');
    elseif ~ischar(odefile) && ~isa(odefile, 'inline')
        error('First argument must be a single-quoted string.  See gni_fgauss1.');
    end

% Test that 'options' is a structure.
    if nargin == 1
        tspan = []; y0 = []; options = [];
    elseif nargin == 2
        y0 = []; options = [];
    elseif nargin == 3
        options = [];
    elseif ~isempty(options) && ~isa(options,'struct')
        error('Correct syntax is gni_gauss(''odefile'',tspan,y0,options).');
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
            msg = sprintf('Use gni_fgauss1(''%s'',tspan,y0,...) instead.',odefile);
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
    
% Get parameters.
    canvector = strcmp(gni_get(options,'Vectorized','off'),'on');
    maxiter = gni_get(options,'MaxIter',50);
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
    if isempty(varargin)
        args = {};
    else
        args = [{[]} varargin];
    end
    if size(y0,2) ~= 1
        y0 = y0';
    end
    ny = length(y0);
    F0 = feval(odefile,t0,y0,args{:});
    nfevals = nfevals + 1;
    [mf,nf] = size(F0);
    if nf > 1
        error([upper(odefile) ' must return a column vector.'])
    elseif mf ~= ny
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
    Ynew = y0;
    if nargout > 0
        nout = 1;
        
        Tout = t0:outstep*h:tf;
        Yout = zeros(ny,length(Tout));
        Yout(:,1) = y0; 
    end
    t = t0;
    
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

% Initialize the IRK method.
    method = gni_get(options,'Method','G12');
    [ns, A, B, C, m, Mu, Alpha, Nu] = coeff_gauss(method);
    uround = 2.221e-16;
    A = h * A';
    B = h * B;
    C = h * C';
    Alpha = h * Alpha;
    Mu = h * Mu';
    Nu = h * Nu';  
    
% preallocation.
    F = zeros(ny,ns);
    Ey = zeros(ny,1);
    ONE = ones(1,ns);
    nsy = ns * ny;

% Simple starting approximation for the first step.
    Zold = F0 * C;

%% THE MAIN LOOP.
    outpoint = 0;
for n = 1:nt-1  
    Yold = Ynew;
    addpoint = false;
    outpoint = outpoint + 1;
    if (outpoint >= outstep)
        outpoint = 0;
        addpoint = true;
    end
    
    dynold = 0;
    dyno = 1;
    niter = 0;
    
    TT = t * ONE + C;
    Yone = Yold * ONE;
    Const = max(0.1,abs(Yone));
% Solve nonlinear system.
    while dyno>=uround    
        %% CALL RK.
        YY = Yone + Zold;
        if (canvector) 
            F = feval(odefile,TT,YY,args{:});
        else            
            for i = 1:ns
                F(:,i) = feval(odefile,TT(i),YY(:,i),args{:});
            end            
        end
        nfevals = nfevals + ns;
        
        Znew = F * A;
        dyno = sqrt(sum(sum(((Zold - Znew)./Const).^2)) / nsy);
        Zold = Znew;   
        
        %% CONVERGENCE TEST.
        niter = niter + 1;
        if (dynold < dyno) && (dyno < 10*uround)
            break;
        end
        if (niter > maxiter)
            if (dyno < 1e5*uround)
                warning('Convergence of fixed-point iteration failed after %d iterations.\nObtained error = %0.5g, continuing anyway...',maxiter,dyno);
                break;          
            else
				error('Convergence of fixed-point iteration failed after %d iterations.\nObtained error = %0.5g, stopping here.',maxiter,dyno);
            end
        end    
        dynold = dyno;             
    end
    niters = niters + niter;
    %% SYSTEM SOLVED.
    
%% LOOP FOR ADVANCING ONE STEP. 
    Ey = Ey + F * B;
    Ynew = Yold + Ey;
    Ey = Ey + (Yold - Ynew);
    
    t = t0 + n*h;
    
% Starting approximation for next step.
    F1 = feval(odefile,t+h,Ynew,args{:});
    Y2 = Yold + [F F0 F1] * Mu;
    F2 = feval(odefile,t+h*m(3),Y2,args{:});
    Zold = F * Alpha + [F0 F1 F2] * Nu;
    F0 = F1;
    
    nfevals = nfevals + 2; 
    
% If we are generating output.
    if addpoint
        if haveoutput
            nout = nout + 1;
            Yout(:,nout) = Ynew;
        end
        if haveoutfun
            if feval(outfun,t,Ynew(outputs),'') == 1
                return;
            end
        end
    end
end
if haveoutfun
  feval(outfun,[],[],'done');
end

%% Problem solved.
% Finalize the output.
    solver_output = {};
    if nargout == 1
        av_niter = niters / n;
        sol.t = Tout;
        sol.y = Yout;
        sol.stats.nfevals = nfevals;
        sol.stats.aviters = av_niter;
        sol.stats.nsteps = n + 1;
        solver_output{1} = sol;
    elseif nargout >= 2
        solver_output{1} = Tout;
        solver_output{2} = Yout;
    end

    varargout = solver_output;
    