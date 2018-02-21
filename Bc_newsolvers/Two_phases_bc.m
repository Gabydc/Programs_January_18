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
for cp = [0]
for per = [6 ]
    for optls = 1
        
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
        
        clearvars -except use_ICCG use_DICCG use_POD optls per cp
        close all
        clc
        mrstModule add incomp
        % Create the grid
        Options
        
        %Initial cell
        [nxi, nyi, nzi] = deal(1, 1, 1);
        if(model_SPE)    % SPE 10
            layers = 1 : 2;
            % Create the grid
            [nx, ny, nz] = deal(60, 220, numel(layers));
        else              % Layered model
            % Create the grid
            [nx, ny, nz] = deal(5, 5,1);
            % Contrast in permeability layers
            %per = 6;
            % Number of layers
            szl = 5;
            % Direction of the layers
            l_dir = 'x';
        end
        N = nx*ny*nz;
        [Lx, Ly, Lz] = deal(nx, ny, nz);
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
        
        gravity on
        maxIterations=500;
        % We define the maximum of iterations for the liner solver
        maxIterations=500;
        % We define the tolerance of the linear solver
        k=8;
        tol = 5*10^(-k);
        tol_p = tol;
        % We can change the permeability one layer is 10*milli*darcy(), the second
        % is 10*milli*darcy()*10^(-per)
        vx = linspace(0, 1, 11) .';
        vy = linspace(1, 0, 11) .';
        pc_form = 'nonwetting';
        [muw, muo] = deal(1,10);
        
        [rhow, rhoo] = deal(1000,700);
        [krw, kro] = deal(2,2);
        cap_scale = 10;
        props = constantProperties([   muw,  muo] .* centi*poise, ...
            [rhow, rhoo] .* kilogram/meter^3);
        [kr, pc]  = tabulatedSatFunc([vx, vx.^krw, vy.^kro, vy.*cap_scale*barsa]);
        % Number of deflation vectors
        dv=10;
        podv{1}=[1:10];
        podv{2}=[6:10];
    
         
                             %% Fluid model
                    % We set up a two-phase fluid. Viscosity, density is
                    % set for oil and water using the valuers from vars
                    Fluid_props

                    
     
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
    units_Pb = meter^3/day;
    % Pressure of the reservoir
    P0 = 100;
end
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
x         = initState(G, [], P0*barsa, [0 1]);
x.wellSol = initWellSol(W, 0);
else
    %x         = initResSol (G, P0*barsa);
    x          = initState(G, [], P0*barsa, [0 1]);
end
%% Solver variables
tol = 8;
tol_p     = 5.0*10^(-tol);
tol_t   = 5.0e-7;
maxIter = 1000;


lsolver = 3;
p0 = x.pressure;
psolverop_1
fn = @(A, b) solver.solveLinearSystem(A, b);
psolve = @(state) incompTPFA_Def(state, G, hT, fluid,  'bc', bc,'MatrixOutput',true,'LinSolve', fn);
[x,preport]    = psolve(x);

plot_props
%% Tempora lparameters
T      = 4800*day;
steps = 240;
DT     = T/steps;

                    
                    %% Deflation parameters
                    % Number of deflation vectors for the window
                    dv    = 10;
                    dp    = 5;
                    if(window)
                        if use_POD
                            dpod   = [dv-dp+1 : dv];
                        else
                            dpod   = [1 : dv];
                        end
                    else    %Training run
                        % POD vectors
                        dpod   = [nstep-dv+1:nstep];
                    end
                    % Last experiment, to close the table
                    last       = false;

    end
end
end