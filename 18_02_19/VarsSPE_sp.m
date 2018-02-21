%% Simulate the SPE10 base case for a single phase
% The simulation uses a DICCG/ ICCG/ agmg pressure solver 

%% Define Solver
use_ICCG   = false;
use_DICCG  = true;
training   = false;
use_agmg   = false;
use_POD    = true;
plot_sol   = true;
save_res   = false;
use_wells  = true;
model_SPE  = true;
use_cp     = false;
window     = true;
% Solvers variables
tol = 5.0e-8;
maxIter = 1500;
dv         = 4;     % Number of deflation vectors for the window
last       = false;  % Last experiment, to close the table
per = 1;
if(window)
    % Deflation parameters
    podv   = 4;                   % Number of POD vectors
    dpod   = [dv-podv+1 : dv];    % POD vectors
else
    %Training run
    dpod   = [nstep-dv+1:nstep];  % POD vectors
end



%% Fluid model
% When gravity forces are absent, the only fluid property we need in the
% incompressible, single-phase flow equation is the viscosity. However, the
% flow solver is written for general incompressible flow and requires the
% evaluation of a fluid object that can be expanded to represent more
% advanced fluid models. Here, however, we only use a simple fluid object
% that requires a viscosity and a density (the latter is needed when gravity
% is present)3 bar
gravity reset off
fluid = initSingleFluid('mu' , 1*centi*poise, ...
                        'rho', 1000*kilogram/meter^3);
%Fluid_props
                    



%% Model
% SPE10 layers
layers = 1 : 1;
% Contrast in permeability layers
%per = 2;

Create_rock; 

%% Compute half transmissibilities
% All we need to know to develop the spatial discretization is the reservoir
% geometry and the petrophysical properties. This means that we can compute
% the half transmissibilities without knowing any details about the fluid
% properties and the boundary conditions and/or sources/sinks that will
% drive the global flow:
T = simpleComputeTrans(G, rock);





%% Initial values
if(use_wells)
    units_P = barsa;
    units_I = barsa;
    %units_I = stb/day;
    %Pressure of the wells
    P = 275;
    P_s = 100;
    I = 1100;
    P0 = 500;
else
    P_b = 40;
    units_Pb = meter^3/day;
    P0 = 100;
end
    % Number of time steps with same pressure (variation in wells)
    tch =2;
    
%% Compute sources
if(use_wells)
    % Set Comp_i=[0,0] in producers to counter X-flow effects...
    
    well_ip = 'ip_tpf';
    W = verticalWell([], G, rock,  1,   1, [], 'Type', 'bhp', ...
        'InnerProduct', well_ip, ...
        'Val', P*barsa, 'Radius', 0.125*meter, ...
        'Name', 'P1');
    
    W = verticalWell(W , G, rock, nx,   1, [], 'Type', 'bhp', ...
        'InnerProduct', well_ip, ...
        'Val', P*barsa, 'Radius', 0.125*meter, ...
        'Name', 'P2');
    
    W = verticalWell(W , G, rock, nx, ny, [], 'Type', 'bhp', ...
        'InnerProduct', well_ip, ...
        'Val', P*barsa, 'Radius', 0.125*meter, ...
        'Name', 'P3');
    
    W = verticalWell(W , G, rock,  1, ny, [], 'Type', 'bhp', ...
        'InnerProduct', well_ip, ...
        'Val', P*barsa, 'Radius', 0.125*meter, ...
        'Name', 'P4');
    
    W = verticalWell(W , G, rock, ceil(nx/2), ceil(ny/2), [], 'Type', 'bhp',   ...
        'InnerProduct', well_ip, ...
        'Val', I*units_I, 'Radius', 0.125*meter, ...
        'Name', 'I1');
    
    %% Changing wells parameters
    

    Changing_w_1ph
else
    %Injection through boundary
    %Boundary conditions. We set boundary conditions,
    %waterflooding, right boundary / wells
    pv = poreVolume(G, rock);
    injRateP = -sum(pv)/(DT*nstep)
    injRate = -(P_b/10)*units_Pb
    bc = fluxside([], G, 'xmin', -injRate, 'sat', [1, 0]);
    bc = pside(bc, G, 'xmax', 0*barsa, 'sat', [0, 1]);
    
end
%%
if(use_wells)
%x         = initResSol (G, P0*barsa);
%x         = initState(G, P0*barsa);

x = initState(G, W, P0*barsa);
x.wellSol = initWellSol(G, 0);
else
% x = initState(G, P0*barsa);
% x         = initResSol (G, P0*barsa);
% %x          = initState(G, P0*barsa);
end
%  state         = initResSol (G, 10);
%    state.wellSol = initWellSol(G, 10);


%% Create directories to save results
%dir = '/mnt/sda2/cortes/Results/2017/Report/12_17/t/';
%dir = '/mnt/sda2/cortes/Results/2017/Report/12_19/check/';
dir = '../2017_CM/Results/';
if(use_wells)
    if(model_SPE)
        folder=[ 'SPE10_' num2str(numel(layers))  'P_' num2str(P)];
    else
        folder=[ 'per_' num2str(per) 'sz_' num2str(nx)  'layers_' num2str(nz)  'P_' num2str(P)];
    end
else
    if(model_SPE)
        folder=[ 'SPE10_' num2str(numel(layers))  'bc_' num2str(P_b)];
    else
        folder=[ 'per_' num2str(per) 'sz_' num2str(nx) '_layers_' num2str(nz)  'bc_' num2str(P_b)];
    end
end
mkdir([dir], folder)
dir1 = [dir folder '/'];

if use_ICCG
    folder=['ICCG' ];
    if training
        folder=['ICCG_t' ];
    end
else
    if(use_wells)
    folder=['DICCG_dv_' num2str(dv) 'pod_' num2str(numel(dpod)) 'Pp_' num2str(P_s)];
    else
    folder=['DICCG_dv_' num2str(dv) 'pod_' num2str(numel(dpod)) 'Pb_' num2str(P_b)];    
    end
end
mkdir([dir1], folder)
dir2 = [dir1 folder '/'];
