% programmer: Gabriela Diaz
% e-mail    : g.b.diazcortes@gmail.com
% date      : 03-02-2018
%% Two-phase reservoir simulation
% This program simulates two-phase flow in heterogeneous porous media.
% The models are:
% The SPE 10 benchmark
% Academic layered problem
% Sequential schemes are used to model the problem, the Saturation is solve
% with MRST, and the solution of the pressure can be obtained with diverse
% methods ('\',agmg,ICCG,DICCG)
% We solve
% $$\nabla\cdot v = q, \qquad
%    v=\textbf{--}\frac{K}{\mu} \bigl[\nabla p+\rho g\nabla z\bigr],$$
% The capillary pressure can be taked into account (cp =1)
for use_Cp = 0
    % The contrast between permeability layers can be varied one layer has
    % permeability of 10*milli*darcy, the secon is varied as
    % 10*10^(-per)*milli*darcy
    for per = [6 ]
        % Choose the solver
        for optls = 2
            
            switch optls
                case 1
                    use_ICCG   = true;
                    use_DICCG  = false;
                    use_POD    = false;
                case 2
                    use_ICCG   = false;
                    use_DICCG  = true;
                    use_POD    = false;
                case 3
                    use_ICCG   = false;
                    use_DICCG  = true;
                    use_POD    = true;
            end
            
            clearvars -except use_ICCG use_DICCG use_POD optls per use_Cp
            close all
            
            
            
            
            
            dpod=[];
            if (use_DICCG)  && (use_POD)
                pod=0;
                continue
            end
            
            %We choose the number of POD deflation vectors to use. This number can be
            %varied in the vars subroutine
            % for rr=1:2
            rr=1
            %                 if (~use_POD) && (rr==2)
            %                     continue
            %                 end
            %                 if (~use_POD) && (use_DICCG)
            %                     continue
            %                 end
            %
            %% Define the model
            % To set up a model, we need: a grid, rock properties (permeability), a
            % fluid object with density and viscosity, and boundary conditions.
            Options
            
            %Initial cell
            [nxi, nyi, nzi] = deal(1, 1, 1);
            if(model_SPE)    % SPE 10
                layers = 1 : 2;
                % Create the grid
                [nx, ny, nz] = deal(60, 220, numel(layers));
            else              % Layered model
                % Create the grid
                [nx, ny, nz] = deal(35, 35,1);
                % Contrast in permeability layers
                %per = 6;
                % Number of layers
                szl = 7;
                % Direction of the layers
                l_dir = 'x';
            end
            N = nx*ny*nz;
            [Lx, Ly, Lz] = deal(35, 35, 1);
            if(model_SPE)
                [cartDims, physDims, G, rock] = Create_rock_f(model_SPE, nx, ny, nz, ...
                    Lx, Ly, Lz, 'layers' , layers);
            else
                [cartDims, physDims, G, rock] = Create_rock_f(model_SPE, nx, ny, nz, ...
                    Lx, Ly, Lz, 'per', per, ...
                    'szl', szl, 'l_dir', l_dir);
                
            end
            
            perm = rock.perm(:,1);
            permD = log10(perm/darcy);
            %% Fluids properties
            Fluid_props
            
            
            %% Temporal variables
            % Time steps
            DT    = 20*day();
            nstep = 100;
            T     = nstep * DT;
            %% Solver variables
            tol = 8;
            tol_v = tol;
            tol_p   = 5.0*10^(-tol);
            tol_t   = 5.0e-7;
            % We define the maximum of iterations for the liner solver
            maxIter = 1000;
            
            %% Deflation parameters
            % Number of deflation vectors for the window
            dv    = 10;
            dp    = 5;
            if(window)
                if use_POD
                    dpod   = [dv-dp+1 : dv];
                    podv{1}=[1:dv];
                    podv{2}=dpod;
                else
                    dpod   = [1 : dv];
                end
            else    %Training run
                % POD vectors
                dpod   = [nstep-dv+1:nstep];
            end
            % Last experiment, to close the table
            last       = false;
            
            
            
            %% Initial values
            % Units of the pressure
            pu = {'barsa', 'stb/day', 'meter^3/day'};
            if(use_wells)
                units_P = pu{1};
                units_I = pu{1};
                % Pressure of the training run
                P = 275;
                if (training)
                    % Pressure of the example
                    P = 400;
                end
                % Pressure of the injector
                I = 1100;
                % Pressure of the reservoir
                P0 = 500;
            else
                %Pressure of the boundary
                P_b = 4;
                % Pressure of the reservoir
                P0 = 100;
            end
            % Number of time steps with same pressure (variation in wells)
            tch =2;
            
            %% Compute sources
            if(use_wells)
                % Include wells
                wtype    = {'bhp', 'bhp', 'bhp', 'bhp', 'bhp'};
                wtarget  = [P * barsa, P * barsa, P * barsa, P * barsa, I * barsa];
                wrad     = [0.125, 0.125, 0.125, 0.125, 0.125] .* meter;
                wloc     = [  nxi,   nxi,     nx,   nx, round(nx/2);
                    nyi,   ny,     nyi,   ny, round(ny/2)];
                wname    = {'P1', 'P2', 'P3', 'P4', 'I'};
                sgn      = [ 1 ,  1 ,  1 ,  1 , -1 ];
                W = [];
                for w = 1 : numel(wtype),
                    W = verticalWell(W, G, rock, wloc(1,w), wloc(2,w), [], ...
                        'Type', wtype{w}, 'Val', wtarget(w), ...
                        'Radius', wrad(w), 'Name', wname{w}, ...
                        'Sign', sgn(w), 'InnerProduct', 'ip_tpf');
                end
                
                % Plot_1(permD, G, 'dnames', {{'Premeability'}},'o_daspect',[6 6 1], ...
                %     'o_ax', 1, 'l_cb', true, 'o_cb', 0, 'W',W, 'figure', 3 ,'o_save', ...
                %     dir,'o_view',[0 90])
                Changing_w
                
            else
                %Injection through boundary
                %Boundary conditions. We set boundary conditions,
                %waterflooding, right boundary / wells
                pv = poreVolume(G, rock);
                %injRate = -sum(pv)/(500*day);
                injRate = -0.4*meter^3/day;
                bc = fluxside([], G, 'xmin', -injRate, 'sat', [1, 0]);
                bc = pside(bc, G, 'xmax', 0*barsa, 'sat', [0, 1]);
                
            end
            
            %% Assemble and solve the linear system
            % To solve the flow problem, we use the standard two-point
            % flux-approximation method (TPFA), which for a Cartesian grid is the same
            % as a classical seven-point finite-difference scheme for Poisson's
            % equation. This is done in two steps: first we compute the
            % transmissibilities and then we assemble and solve the corresponding
            % discrete system.
            hT   = computeTrans(G, rock);
            
            % Initialize the reservoir
            if(use_wells)
                %x         = initResSol (G, P0*barsa);
                rSol         = initState(G, [], P0*barsa, [0 1]);
                rSol.wellSol = initWellSol(W, 0);
            else
                %x         = initResSol (G, P0*barsa);
                rSol          = initState(G, [], P0*barsa, [0 1]);
            end
            
            IWS =  rSol.s(:,1);
            IWP = rSol.pressure;
            % Set linear solver for incompTPFA_Def ( The difference with incompTPFA is
            % that we can get the report
            p0 = rSol.pressure;
            %Select a linear solver 1.AGMG, 2.GMRES, 3.PCG_IC, 4.DPCG_IC, other Backslash
            lsolver = 3;
            % Select transport solver 1. Explicit, 2. Implicit
            tsolver=2;
            psolverop_1
            solver = ICCGSolverAD('tolerance', tol_p,'maxIterations',  maxIter,'x0',p0);
            fn = @(A, b) solver.solveLinearSystem(A, b);
            psolve = @(state) incompTPFA_Def(state, G, hT, fluid, 'bc', bc ,'MatrixOutput',true,'LinSolve', fn);
            
            
            %% Set Transport Solver
            
            linsolve_t = @(J, F) agmg(J, F, 50, tol_t, maxIter, 0);
            
            
            
            
            %% Define transport solver
            switch tsolver
                case 1
                    % Explicit tranport solver
                    if(use_wells)
                        tsolve = @(x, dt, fluid) ...
                            explicitTransport(x, G, dt, rock, fluid, 'wells', W, ...
                            'LinSolve', linsolve_t, 'verbose',false);
                    else
                        tsolve = @(x, dt, fluid) ...
                            explicitTransport(x, G, dt, rock, fluid, 'bc', bc, ...
                            'LinSolve', linsolve_t, 'verbose',false);
                    end
                    % tsolve  = @(state, dT, fluid) ...
                    %     explicitTransport(state, G, dT, rock, fluid,  'bc', bc,  'verbose', false);
                case 2
                    % Implicit transport solver: try with one time step
                    if(use_wells)
                        tsolve = @(x, dt, fluid) ...
                            implicitTransport(x, G, dt, rock, fluid, 'wells', W, ...
                            'LinSolve', linsolve_t, 'verbose',false);
                    else
                        tsolve = @(x, dt, fluid) ...
                            implicitTransport(x, G, dt, rock, fluid, 'bc', bc, ...
                            'LinSolve', linsolve_t, 'verbose',false);
                    end
                    %  tsolve  = @(state, dT, fluid) implicitTransport(state, G, dT, rock, fluid,  'bc', bc, 'Verbose', false);
            end
            
            
            
            
            %% Initiate pressure solver
            % Solve initial pressure in reservoir
            rSol = psolve(rSol);
            A01 = rSol.A;
            
            
            %Create the directory
            create_dir
            if (use_DICCG) && (~window)
                use_DICCG
                filep=['Pressure'];
                filename=[dir1 filep];
                load(filename)
                [U,S]=PODbasis(Pressure);
                Z = U(:,dpod);
            end
            podi = 0;
            
            
            % If we use wells, we can define here the position of the wells
            % w1=10*2^pw;
            % w2=2*w1;
            % Set up grid and petrophysical data (varsbc)
            % We use a Cartesian grid of size nx-by-ny with homogeneous petrophysical
            % data: permeability of 100 mD and porosity of 0.2.
            
            dTplot = ceil(T/3);  % plot only every T/3 day
            ts=0;
            podi=0;
            t =0;
            % while t < T,
            for k=1 : 11;
                
                
                k
                
                ts=ts+1;
                
                % TRANSPORT SOLVER
                t0 = tic;
                [rSol,treport(ts)]   = tsolve(rSol, DT, fluid);
                
                dt_t(k) = toc(t0);
                x_s1 = rSol.s;
                save x_s1 x_s1
                load x_s
                %                             comparevars1(x_s,x_s1,'Saturatio_1')
                %                     break
                % Check for inconsistent saturations
                s = [rSol.s(:,1)];
                assert(max(s) < 1+eps && min(s) > -eps);
                % PRESSURE SOLVER
                p0 = rSol.pressure;
                if(use_wells)
                    for i=1:numel(W)
                        p0(N+i) = 0;
                    end
                    if(training)
                        W = W1{k};
                        %  W(5).val = I_P(k);
                    else
                        W = W0{k};
                        %   W(5).val = I_P(k);
                    end
                end
                if  use_ICCG
                    use_ICCG
                    
                    if(use_wells)
                        if training
                            W = W1{k};
                            %  W(5).val = I_P(k);
                        else
                            W = W0{k};
                            %   W(5).val = I_P(k);
                        end
                        %   W(5).val = PI_1;
                        solver = ICCGSolverAD('tolerance', tol_p,'maxIterations',  maxIter,'x0',p0,'W', W);
                        linsolve_p = @(A, b) solver.solveLinearSystem(A, b);
                        psolve = @(x) incompTPFA_Def(x, G, hT, fluid, 'wells', W,'LinSolve', linsolve_p);
                    else
                        [psolve]= ICCG_b(tol_p,maxIter,p0,G,hT,fluid,bc,true);
                    end
                    
                    t0 = tic;
                    [rSol,preport(k)]= psolve(rSol);
                    dt_p(k) = toc(t0);
                    
                    %
                elseif (use_DICCG)
                    use_DICCG
                    if (window)
                        if k < dv+1
                            use_ICCG
                            if(use_wells)
                                solver = ICCGSolverAD('tolerance', tol_p,'maxIterations',  maxIter,'x0',p0,'W', W);
                                linsolve_p = @(A, b) solver.solveLinearSystem(A, b);
                                psolve = @(x) incompTPFA_Def(x, G, hT, fluid, 'wells', W,'LinSolve', linsolve_p);
                            else
                                [psolve]= ICCG_b(tol_p,maxIter,p0,G,hT,fluid,bc,true);
                            end
                            t0 = tic;
                            [rSol,preport(k)]= psolve(rSol);
                            dt_p(k) = toc(t0);
                            POD_V(:,k) = rSol.pressure;
                            
                        else
                            %                 % Number of POD vectors
                            %                 if use_POD
                            %                     dpod   = [k-dp : k-1];
                            %                 else
                            %                     dpod   = [k-dv : k-1];
                            %                 end
                            
                            %                 podi=podi+1;
                            %                 Z=POD_V(:,podi:podi+dv-1);
                            %                 if use_POD
                            %
                            %                     [U,S]=defpodf_Dt(Z,dir2,dv,t/day(),dTplot/day());
                            %                     Z=U(:,dpod);
                            %                 end
                            
                            [U,S]=PODbasis(POD_V(:,k-dv:k-1));
                            Z = U(:,dpod);
                            if(use_wells)
                                for i=1:numel(W)
                                    Z(N+i,1)=0;
                                end
                                
                            end
                            if(use_wells)
                                W = W0{k};
                                %  W(5).val = I_P(k);
                                solver = DICCGSolverAD('tolerance', tol_p,'maxIterations',  maxIter,'Z',Z,'x0',p0,'W', W,'dir',dir2);
                                linsolve_p = @(A, b) solver.solveLinearSystem(A, b);
                                psolve = @(x) incompTPFA_Def(x, G, hT, fluid, 'wells', W,'LinSolve', linsolve_p,'MatrixOutput',true);
                            else
                                solver = DICCGSolverAD('tolerance', tol_p,'maxIterations',  maxIter,'Z',Z,'x0',p0,'bc', bc,'dir',dir2);
                                linsolve_p = @(A, b) solver.solveLinearSystem(A, b);
                                psolve = @(x) incompTPFA_Def(x, G, hT, fluid, 'bc', bc,'LinSolve', linsolve_p,'MatrixOutput',true);
                            end
                        end
                    else
                        if(use_wells)
                            W = W0{k};
                            for i=1:numel(W)
                                Z(N+i,1)=0;
                            end
                            %  W(5).val = I_P(k);
                            solver = DICCGSolverAD('tolerance', tol_p,'maxIterations',  maxIter,'Z',Z,'x0',p0,'W', W);
                            linsolve_p = @(A, b) solver.solveLinearSystem(A, b);
                            psolve = @(x) incompTPFA_Def(x, G, hT, fluid, 'wells', W,'LinSolve', linsolve_p);
                        else
                            solver = DICCGSolverAD('tolerance', tol_p,'maxIterations',  maxIter,'Z',Z,'x0',p0,'bc', bc);
                            linsolve_p = @(A, b) solver.solveLinearSystem(A, b);
                            psolve = @(x) incompTPFA_Def(x, G, hT, fluid, 'bc', bc,'LinSolve', linsolve_p);
                        end
                    end
                    t0 = tic;
                    [rSol,preport(k)]= psolve(rSol);
                    dt_p(k) = toc(t0);
                    
                    if(window)
                        POD_V(:,k) = rSol.pressure;
                    end
                end
                %                         Zp(:,ts)=rSol.pressure;
                %                         if ts <dv+1
                %                             % Update solution of pressure equation. If backslash is used, there
                %                             % is no report
                %                             [rSol,preport(ts)]    = psolve(rSol);
                %                         else
                %                             podi=podi+1;
                %                             Z=Zp(:,podi:podi+dv-1);
                %                             if use_POD
                %
                %                                 [U,S]=defpodf_Dt(Z,dir2,dv,t/day(),dTplot/day());
                %                                 Z=U(:,dpod);
                %                             end
                %                             lsolver = 4;
                %                             p0 = rSol.pressure;
                %                             psolverop_1
                %                             fn = @(A, b) solver.solveLinearSystem(A, b);
                %                             psolve = @(state) incompTPFA_Def(state, G, hT, fluid,  'bc', bc,'MatrixOutput',true,'LinSolve', fn);
                %                             [rSol,preport(ts)]    = psolve(rSol);
                %
                %                         end
            end
            t = t + DT;
            Sat_w(:,ts) =  rSol.s(:,1);
            Pressure1(:,ts) = rSol.pressure(:)/barsa;
            
            % end
            
            
            
            
        end
    end
end
def =use_DICCG;
pod = use_POD;
savefilesf
plot_props

xp1 = rSol.pressure;
%  save xs1 xs1
load xp
comparevars1(xp,xp1,'Pressure')
x_s1 = rSol.s;
% save x_s1 x_s1
load x_s
comparevars1(x_s,x_s1,'Saturation')
load A0
comparevars(A0,A01,'A0')

A1 = rSol.A;
a=1
load A
comparevars(A,A1,'A')
load U_1
comparevars(U,U_1,'eigvec')
load S1
comparevars1(S,S1','eigval')
%%
for i = 1 :10
    figure
    plot(U(:,i),'or')
    hold on
    plot(U_1(:,i),'*b')
    pause
end

% A = rSol.A;
% save A A
% clear all
% close all
% load A0
% load A01
% comparevars(A0,A01)
% load A
% load A1
% comparevars(A,A1)
