%...  The MatMol Group (2016)
    function varargout = gni_flmm2(odefile,tspan,y0,options,varargin)


%GNI_FLMM2  Solve non-stiff differential equations.
%   [T,Q,P] = GNI_FLMM2('G',TSPAN,[Q0;P0]) with TSPAN = [T0 TFINAL] 
%   integrates the system of differential equations q'' = g(t,q) 
%   from time T0 to TFINAL with initial conditions q = Q0, q' = P0. 'G' is
%   a scalar. For a vector Y, G(~,Y) must return a column vector
%   corresponding to g(y). Each column in the solution array Y corresponds 
%   to a time returned in the column vector T.
%   
%   Use [T,Q] = GNI_FLMM2('G',TSPAN,[Q0;P0],OPTIONS) if you don't use P.
%
%   [T,Q,P] = GNI_FLMM2('G',TSPAN,[Q0;P0],OPTIONS) solves as above 
%   with default integration properties replaced by values in OPTIONS, 
%   an argument created with the GNI_SET function. Commonly used options 
%   are step size 'StepSize' (1e-2 by default) and method used 'Method'
%   ('8b' by default).  
%
%   SOL = GNI_FLMM2(ODEFILE,TSPAN,Y0,...) returns a structure. The
%   structure SOL includes these fields:
%   sol.t          Steps chosen by the solver
%   sol.y          Each column sol.y(:,i) contains the solution at sol.t(i)
%   sol.solver     Solver name  
%   sol.stats      Statistics
%
%   Example   
%       sol = gni_flmm2('G_solar',[0 1000000]);   
%       plot_solar(sol.q)

%     solves the solar system q'' = G_solar(t,y), over 1 billon years using 
%     initial values and step size defined in the file 'G_solar' and the
%     default method '8b' and plots the components of the solution in the
%     phase space with the help of the function 'plot_solar'.
%
%   gni_flmm2 is an implementation of fixed symmetric multistep methods. 
%   It is for autonomous second order Hamiltonian systems.
%   
%   REFERENCE:
%   E. Hairer, C. Lubich, G. Wanner, Geometric Numerical Integration. 
%   Structure-Preserving Algorithms for Ordinary Di?erential Equations,
%   springer-verlag, berlin, second edition Edition, Vol. 31, Springer
%   Series in Computational Mathematics 31, 2006.

    solver_name = 'gni_flmm2';

% stats
    nfevals = 0; 

% Test that 'odefile' is a string
    if nargin == 0
        error('Not enough input arguments.  See gni_flmm2.');
    elseif ~ischar(odefile) && ~isa(odefile, 'inline')
        error('First argument must be a single-quoted string.  See gni_flmm2.');
    end
    
% Test that 'options' is a structure
    if nargin == 1
        tspan = []; y0 = []; options = [];
    elseif nargin == 2
        y0 = []; options = [];
    elseif nargin == 3
        options = [];
    elseif ~isempty(options) && ~isa(options,'struct')
        error('Correct syntax is gni_flmm2(''odefile'',tspan,y0,options).');
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
    tf = tspan(2); 
    if t0 >= tf
        error('The final time must greater than the starting time.');
    end
    t = t0;
    
% get parameters
    canvector = strcmp(gni_get(options,'Vectorized','off'),'on');
    compensated = strcmp(gni_get(options,'Compensated','on'),'on');
    maxiter = gni_get(options,'MaxIter',50);
    h = gni_get(options,'StepSize');
    hold = h;
    if (isempty(h))
        nsteps = gni_get(options,'NumSteps');
        if (isempty(nsteps))
            h = 1e-2;
            fprintf('Warning: No initial step size provided, using h=1e-2 instead.\n');
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
    
% Initialise the method
    method = gni_get(options,'Method','8b');
    [k,A,B,C] = coeff_lmm(method);
    A = A';
    hB = h * B';
    Ch = C' / h;
    
% Test that the length of y0 is pair
    if size(y0,2) ~= 1
        y0 = y0';
    end
    ny = length(y0);
    if (mod(ny,2) ~= 0)
        error('The initial value must have 2*N components, where N is the dimension of the system.');
    end
    npq = ny / 2;
    
% Test that odefile is consistent.
    Qn = y0(1:npq);
    Pn = y0(npq+1:ny);
    if isempty(varargin)
        args = {};
    else
        args = [{[]} varargin];
    end
    FS = feval(odefile,[],Qn,args{:}); 
    nfevals = nfevals + 1;
    [mf,nf] = size(FS);
    if nf > 1
        error([upper(odefile) ' must return a column vector.'])
    elseif mf ~= npq
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

%% Compute first step by a Gauss-legendre method
% Initialise the IRK
    uround = 2.221e-16;    
    [ns,Ar,Br,Cr,nm,mu,alpha,nu] = coeff_gauss('G14');    
    Ar = Ar';
    AA = h^2 * (Ar*Ar);
    BB = h^2 * (Ar*Br);
    Br = h * Br;
    C = h * Cr';
    Alpha = h * alpha;
    Mu = h * mu';
    Nu = h * nu';
% Preallocation   
    KK = zeros(npq,ns);
    Ep = zeros(npq,1);
    Eq = zeros(npq,1);  
    F = zeros(npq,k-1);
    ONE = ones(1,ns);
    QE = zeros(npq,k+1);
    QE(:,1) = Qn;
    Z = zeros(npq,k);
% Starting approximation for the first step
    Zold = Pn * C + 0.5 * FS * C.^2;
    PS = Pn;    
    
% start IRK
    outpoint = 0;
for n=1:k-1
    Qold = Qn;
    Pold = Pn;

    if (n > 1)
        QH = Zold * Mu(1:ns) + Mu(ns+1)*PS + Mu(ns+2)*Pold + Qold;
		Zold = KK * Alpha + FS * Nu(1,:);
		FS = feval(odefile,t+h,Qold,args{:});
		F(:,1) = feval(odefile,t+h*nm(3),QH,args{:});
		PS = Pold;
		Zold = Zold + FS * Nu(2,:) + F(:,1) * Nu(3,:) + Pold * C;	        
        F(:,n)= FS;
        
        nfevals = nfevals + 2;
    end
    
    dynold = 0;
    dyno = 1;
    niter = 0;
    while (dyno > uround)
        QQ = Qold * ONE + Zold;
        if (canvector)
            KK = feval(odefile,t*ONE + C,QQ,args{:});
        else
            for i=1:ns
                KK(:,i) = feval(odefile,t + C(i),QQ(:,i),args{:});
            end
        end
        nfevals = nfevals + ns;
        
        Znew = Pold * C + KK*AA;        
		dyno = sqrt(sum(sum(((Zold - Znew)./max(0.1,abs(Qold)*ONE)).^2)) / (ns*npq));
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
    
    Ep = Ep + KK * Br;
    Pn = Pold + Ep;
    Ep = Ep + (Pold - Pn);
    
    Eq = Eq + h * Pold + KK * BB;
    Qn = Qold + Eq;
    Eq = Eq + (Qold - Qn);

    QE(:,n+1) = Qn;
    Z(:,n+1) = Pn;
end

for position=2:k/2
    addpoint = false;
    outpoint = outpoint + 1;
    if (outpoint >= outstep)
        outpoint = 0;
        addpoint = true;
    end
    if addpoint
        if haveoutput
            nout = nout + 1;
            Qout(:,nout) = QE(:,position);
            Pout(:,nout) = Z(:,position);
        end
        if haveoutfun
            t = t0 + position*h;
            Yn = [QE(:,position);Z(:,position)];
            if feval(outfun,t,Yn(outputs),'') == 1
                return;
            end
        end
    end
end


%% Start lmm
    
% Indices used    
    k2 = k / 2;
    k21 = k2 + 1;
    k1 = k + 1;
    I1_k = 1:k;
    I1_1k = 1:k-1;
    I1_2k = 1:k-2;
    I1_k2 = 1:k2;
    I2_k = 2:k;
    I2_k1 = 2:k+1;
    I2_1k = 2:k-1;
    I2_k2 = 2:k2;
    I1k_k2 = k-1:-1:k2;
    I1k_k21 = k-1:-1:k2+1;
    Ik1_k22 = k+1:-1:k2+2;

% Initialisation
    F(:,I1_2k) = F(:,I2_1k);
    Z(:,I1_1k) = (QE(:,I2_k) - QE(:,I1_1k))/h;
    Ek = zeros(npq,k);
    E = 0;

% THE MAIN LOOP
    for i=k+1:nt+k2
        addpoint = false;
        outpoint = outpoint + 1;
        if (outpoint >= outstep)
            outpoint = 0;
            addpoint = true;
        end        
        
        F(:,k-1) = feval(odefile,[],QE(:,k),args{:});
        if (compensated)
            s1 = (F(:,I1_k2) + F(:,I1k_k2)) * hB;
            s2 = (Z(:,I2_k2) - Z(:,I1k_k21)) * A;
            d = (Ek(:,I2_k2) - Ek(:,I1k_k21)) * A;
            a = Z(:,1);
            Ek(:,k) = s1 + s2 + d + Ek(:,1);
            Z(:,k) = a + Ek(:,k);
            Ek(:,k) = (a - Z(:,k)) + Ek(:,k);
            Ek(:,I1_1k) = Ek(:,I2_k);
            b = QE(:,k);
            E = h * Z(:,k) + E;
            QE(:,k1) = b + E;
            E = (b - QE(:,k1)) + E;
        else
            Z(:,k) = Z(:,1)+(F(:,I1_k2)+F(:,I1k_k2))*hB + ...
                (Z(:,I2_k2)-Z(:,I1k_k21))*A;
            QE(:,k1) = QE(:,k) + h * Z(:,k);
        end
        
    % If we are generating output.
        if addpoint
            Qn = QE(:,k21);
            Pn = (QE(:,Ik1_k22) - QE(:,I1_k2))*Ch;
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
        
        F(:,I1_2k) = F(:,I2_1k);
        Z(:,I1_1k) = Z(:,I2_k);
        QE(:,I1_k) = QE(:,I2_k1);
    end

    nfevals = nfevals + nt - k;
    
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
        sol.stats.nfevals = nfevals;
        sol.stats.nsteps = i;
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
    