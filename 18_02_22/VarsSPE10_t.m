%% Simulate the SPE10 base case
% The simulation uses a DICCG/ ICCG/ agmg pressure solver and an implicit transport
% solver.
%% Time
% Time steps
DT    = 100*day;
nstep =  80;
%% Define Solver
use_ICCG   = true;
use_DICCG  = false;
training   = true;
use_agmg   = false;
use_POD    = true;
plot_sol   = true; 
save_res   = true;
use_wells  = false;

% Solvers variables
tol = 5.0e-8;
maxIter = 1500;

% Deflation parameters
dv = 20;
dpod = [nstep-dv+1:nstep];
last       = false;
%% Fluids properties
fluid = initSimpleFluid('mu' , [   1,  10]*centi*poise     , ...
                        'rho', [1000, 700]*kilogram/meter^3, ...
                        'n'  , [   2,   2]);

%% Initial values
if(use_wells)
    units_P = barsa;
    units_I = barsa;
    %units_I = stb/day;
    %Pressure of the wells
    P = 600;
    I = 1100;
else
    P_b = 10;
    units_Pb = meter^3/day;
end



%% Model

layers = 1 : 1;
[nx, ny, nz] = deal(60, 220, numel(layers));
cartDims = [nx, ny,nz];
physDims = cartDims .* [20, 10, 2]*ft;
N = nx*ny*nz;

% Construct the model
rock     = SPE10_rock(layers);
rock.perm = convertFrom(rock.perm, milli*darcy);
is_pos             = rock.poro > 0;
rock.poro(~is_pos) = min(rock.poro(is_pos));
G = computeGeometry(cartGrid(cartDims, physDims));
T = computeTrans(G, rock);


%% Compute sources
if(use_wells)
    % Set Comp_i=[0,0] in producers to counter X-flow effects...
    
    well_ip = 'ip_tpf';
    W = verticalWell([], G, rock,  1,   1, [], 'Type', 'bhp', ...
        'InnerProduct', well_ip, ...
        'Val', P*barsa, 'Radius', 0.125*meter, ...
        'Name', 'P1', 'Comp_i', [0, 0]);
    
    W = verticalWell(W , G, rock, nx,   1, [], 'Type', 'bhp', ...
        'InnerProduct', well_ip, ...
        'Val', P*barsa, 'Radius', 0.125*meter, ...
        'Name', 'P2', 'Comp_i', [0, 0]);
    
    W = verticalWell(W , G, rock, nx, ny, [], 'Type', 'bhp', ...
        'InnerProduct', well_ip, ...
        'Val', P*barsa, 'Radius', 0.125*meter, ...
        'Name', 'P3', 'Comp_i', [0, 0]);
    
    W = verticalWell(W , G, rock,  1, ny, [], 'Type', 'bhp', ...
        'InnerProduct', well_ip, ...
        'Val', P*barsa, 'Radius', 0.125*meter, ...
        'Name', 'P4', 'Comp_i', [0, 0]);
    
    W = verticalWell(W , G, rock, ceil(nx/2), ceil(ny/2), [], 'Type', 'bhp',   ...
        'InnerProduct', well_ip, ...
        'Val', I*units_I, 'Radius', 0.125*meter, ...
        'Name', 'I1', 'Comp_i', [1, 0]);
    
    %% Changing wells parameters
    
    % Number of time steps with same pressure
    tch =2;
    Changing_w
else
    %Injection through boundary
    %Boundary conditions. We set boundary conditions,
    %waterflooding, right boundary / wells
    pv = poreVolume(G, rock);
    %injRate = -sum(pv)/(500*day);
    injRate = -P_b/10*units_Pb;
    bc = fluxside([], G, 'xmin', -injRate, 'sat', [1, 0]);
    bc = pside(bc, G, 'xmax', 0*barsa, 'sat', [0, 1]);
    
end
%%
if(use_wells)
x         = initResSol (G, 0);
x.wellSol = initWellSol(W, 0);
else
x         = initResSol (G, 0);
end

%%
linsolve_p = @(S, h) agmg(S, h,  1, 5.0e-11, 1000, 0);
linsolve_t = @(J, F) agmg(J, F, 50, 5.0e-11, 2000, 0);

if(use_wells)
   psolve = @(x) ...
      incompTPFA(x, G, T, fluid, 'wells', W, 'LinSolve', linsolve_p);
   tsolve = @(x, dt) ...
      implicitTransport(x, G, dt, rock, fluid, 'wells', W, ...
                        'LinSolve', linsolve_t, 'verbose',true);
else
   psolve = @(x) ...
      incompTPFA(x, G, T, fluid, 'bc', bc, 'LinSolve', linsolve_p);
   tsolve = @(x, dt) ...
      implicitTransport(x, G, dt, rock, fluid, 'bc', bc, ...
                        'LinSolve', linsolve_t, 'verbose',true);    
end




%% Create directories to save results
dir = '/mnt/sda2/cortes/Results/2017/Report/SPE10/training/12_01/ex1/';
%dir = '../SPE10_results/';
if(use_wells)
folder=[ 'SPE10_' num2str(numel(layers))  'DT_' num2str(DT/day) 'step_' num2str(nstep) 'P_' num2str(P)];
else
folder=[ 'SPE10_' num2str(numel(layers))  'DT_' num2str(DT/day) 'step_' num2str(nstep) 'bc_' num2str(P_b)];   
end
mkdir([dir], folder)
dir1 = [dir folder '/'];

if use_ICCG
    folder=['ICCG' ];
    if training
        folder=['ICCG_t' ];
    end
else
    folder=['DICCG' num2str(dv)];
end
mkdir([dir1], folder)
dir2 = [dir1 folder '/'];
