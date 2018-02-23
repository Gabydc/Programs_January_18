%...  The MatMol Group (2016)
    function varargout = gni_vgauss1(odefile,tspan,y0,options,varargin)
%... 
%GNI_VGAUSS1  Solve non-stiff autonomous differential equations.
%   [T,Y] = GNI_VGAUSS1('F',TSPAN,Y0) with TSPAN = [T0 TFINAL] integrates 
%   the autonomous system of differential equations y' = f(y) from time T0 
%   to TFINAL with initial conditions Y0. 'F' is a string. For a scalar T 
%   and a vector Y, F(T,Y) must return a column vector corresponding to 
%   f(y). Each column in the solution array Y corresponds to a time 
%   returned in the column vector T.
%   
%   [T,Y] = GNI_VGAUSS1('F',TSPAN,Y0,OPTIONS) solves as above with default
%   integration properties replaced by values in OPTIONS, an argument 
%   created with the GNI_SET function. Commonly used options are step size 
%   'Precision' (5e-2 by default) and method used 'Method' (G10 by default).  
%
%   Example    
%         [t,y]=gni_vgauss1('F_kepler',[0 20],[0.4 0 0 2]);   
%         plot(y(1,:),y(2,:));
%     solves the system y' = F_kepler(t,y), using the default parameter of
%     precision 5e-2 and the method 'G10' and plots the two first component
%     of the solution. 
%
%   gni_vgauss1 is an implementation of the Gauss methods with variable 
%   step size. It uses fixed point iteration to solve the nonlinear 
%   system at every step. It is for non-stiff Hamiltonian and/or reversible
%   systems.
%   
%   REFERENCE:
%   E. Hairer, C. Lubich, G. Wanner, Geometric Numerical Integration. 
%   Structure-Preserving Algorithms for Ordinary Di?erential Equations,
%   springer-verlag, berlin, second edition Edition, Vol. 31, Springer
%   Series in Computational Mathematics 31, 2006.

solver_name = 'gni_vgauss1';

% Stats.
    nfevals = 0; 
    niters = 0;
    
% Test that 'odefile' is a string.
    if nargin == 0
        error('Not enough input arguments.  See gni_vgauss1.');
    elseif ~ischar(odefile) && ~isa(odefile, 'inline')
        error('First argument must be a single-quoted string.  See gni_vgauss1.');
    end

% Test that 'options' is a structure.
    if nargin == 1
        tspan = []; y0 = []; options = [];
    elseif nargin == 2
        y0 = []; options = [];
    elseif nargin == 3
        options = [];
    elseif ~isempty(options) && ~isa(options,'struct')
        error('Correct syntax is gni_vgauss1(''odefile'',tspan,y0,options).');
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
            msg = sprintf('Use gni_gauss(''%s'',tspan,y0,...) instead.',odefile);
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
    tfin = tspan(2); 
    if t0 >= tspan(2)
        error('The final time must greater than the starting time.');
    end
    t = t0;

% Get parameters.
    canvector = strcmp(gni_get(options,'Vectorized','off'),'on');
    maxiter = gni_get(options,'MaxIter',50);
    epsi = gni_get(options,'Precision',0.05);
    if ((~isa(epsi,'double')) || (size(epsi,1) ~= 1) || (size(epsi,2) ~= 1)) 
        error('The option ''Precision'' must contain a single number');
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
    F0 = feval(odefile,[],y0,args{:});  % if F0=0 ?? 
    nfevals = nfevals + 1;
    [mf,nf] = size(F0);
    if nf > 1
        error([upper(odefile) ' must return a column vector.'])
    elseif mf ~= ny
        msg = sprintf('an initial condition vector of length %d.',mf);
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
    Yn = y0;
    if nargout > 0
        nout = 1;
        
        Yout = zeros(ny,1);
        Yout(:,1) = y0; 
        Tout(:,1) = t;
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
        feval(outfun,[t0 tfin],y0(outputs),'init');
    end

% Initialize the IRK method.
    method = gni_get(options,'Method','G10');
    [ns, A, B, C, ~, mu, alpha, nu] = coeff_gauss(method);
    uround = 2.221e-16;
    hA = epsi * A';
    hB = epsi * B';
    hC = epsi * C';
    hAlpha = epsi * alpha;
    hMu = epsi * mu';
    hNu = epsi * nu';
    
% Preallocation.
    F = zeros(ny,ns);
    sig = zeros(ns,1);
    EY = zeros(ny,1);
    ONE = ones(1,ns);
    nsy = ns * ny;
    
% Simple starting approximation for first step.
    Zold = F0 * hC / norm(F0);
    
%% THE MAIN LOOP.
    outpoint = 0;
    Told = t;
    Yold = Yn;
    done = true;  
    n = 1;
while done  
    addpoint = false;
    outpoint = outpoint + 1;
    if (outpoint >= outstep)
        outpoint = 0;
        addpoint = true;
    end
    
    dynold = 0;
    dyno = 1;
    niter = 0;
    Yone = Yold * ONE;
% Solve nonlinear system.
    Const = max(0.1,abs(Yold)*ONE);
    while dyno>=uround
        %% CALL RK.
        YY = Yone + Zold;
        if (canvector)
            Fsauv = feval(odefile,[],YY,args{:});
            for i = 1:ns
                sig(i,1) = norm(Fsauv(:,i));
                F(:,i) = Fsauv(:,i) / sig(i,1);
            end
        else
            for i = 1:ns
                Fsauv = feval(odefile,[],YY(:,i),args{:});
                sig(i,1) = norm(Fsauv);
                F(:,i) = Fsauv / sig(i,1);
            end
        end
        nfevals = nfevals + ns;
        
        Znew = F * hA;
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
    
    % LOOP FOR ADVANCING ONE STEP.
        EY = EY + F * hB';
        Yn = Yold + EY;
        EY = EY + (Yold - Yn);

        t = Told + hB * (1 ./ sig);
    
    % ADJUST LAST STEP IF OVERFLOW.
        if t >= tfin  
            done = false;
            if t > tfin
                epsi = tfin - Told;
                hA = epsi * A';
                hB = epsi * B';
                hC = epsi * C';
                F0 = feval(odefile,[],Yold,args{:});
                nfevals = nfevals + 1;
                Zold = F0 * hC;
                dynold = 0;
                dyno = 1;
                niter = 0;
                while dyno>=uround
                    YY = Yold * ONE + Zold;
                    if (canvector)
                        F = feval(odefile,[],YY,args{:});
                    else
                        for i = 1:ns
                            F(:,i)= feval(odefile,[],YY(:,i),args{:});
                        end
                    end
                    nfevals = nfevals + ns;
                    Znew = F * hA;
                    dyno = sqrt(sum(sum(((Zold - Znew)./max(0.1,abs(Yold)*ONE)).^2)) / (nsy));
                    Zold = Znew;                 
                    niter = niter + 1;
                    if (dynold < dyno) && (dyno < 10*uround)
                        break;
                    end
                    if (niter > maxiter)
                        if (dyno < 1e5*uround)
                            warning('Convergence of Gauss-Newton failed after %d iterations.\nObtained error = %0.5g, continuing anyway...',maxiter,dyno);
                            break;
                        else
                            error('Convergence of Gauss-Newton failed after %d iterations.\nObtained error = %0.5g, stopping here.',maxiter,dyno);
                        end
                    end
                    dynold = dyno;
                end
                Yn = Yold + F * hB';
                t = tfin;
            end
        else
        % STARTING APPROXIMATION for next step.
            F1 = feval(odefile,[],Yn,args{:});
            F1 = F1 / norm(F1);
            Y2 = Yold + [F F0 F1] * hMu;
            F2 = feval(odefile,[],Y2,args{:});
            F2 = F2 / norm(F2);
            Zold = F * hAlpha + [F0 F1 F2] * hNu;
            F0 = F1;  
            
            nfevals = nfevals + 2; 
        end
        
    % If we generate output.
        if addpoint
            if haveoutput
                nout = nout + 1;
                Tout(nout) = t;
                Yout(:,nout) = Yn;
            end
            if haveoutfun 
                if feval(outfun,t,Yn(outputs),'') == 1
                    return;
                end
            end            
        end
     
     % Advance the integration one step.
        Told = t;
        Yold = Yn;
        
        n = n + 1;      
end
%% Problem solved

if haveoutfun
    feval(outfun,[],[],'done');
end

% Finalize the output.
    solver_output = {};
    if nargout == 1
        av_niter = niters / (n-1);
        sol.t = Tout;
        sol.y = Yout;
        sol.stats.nfevals = nfevals;
        sol.stats.aviters = av_niter;
        sol.stats.nsteps = n;
        solver_output{1} = sol;
    elseif nargout >= 2
        solver_output{1} = Tout;
        solver_output{2} = Yout;
    end

    varargout = solver_output;
    