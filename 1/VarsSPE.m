%% Simulate the SPE10 base case
% The simulation uses a DICCG/ ICCG/ agmg pressure solver and an implicit transport
% solver.
%% Time
% Time steps
DT    = 1*day;
nstep =  150;
%% Define Solver
use_ICCG   = false;
use_DICCG  = true;
training   = false;
use_agmg   = false;
use_POD    = true;
plot_sol   = false;
save_res   = true;
use_wells  = true;
model_SPE  = true;
use_cp     = false;
window     = true;
% Solvers variables
tol = 5.0e-8;
maxIter = 1500;
dv         = 10;     % Number of deflation vectors for the window
last       = false;  % Last experiment, to close the table
per = 1;
if(window)
    % Deflation parameters
    podv   = dv;                   % Number of POD vectors
    dpod   = [dv-podv+1 : dv];    % POD vectors
else
    %Training run
    dpod   = [nstep-dv+1:nstep];  % POD vectors
end



%% Fluids properties
Fluid_props
                    



%% Model
% SPE10 layers
layers = 1 : 85;
% Contrast in permeability layers
%per = 2;

Create_rock; 


T = computeTrans(G, rock);

gravity on


%% Initial values
if(use_wells)
    units_P = barsa;
    units_I = barsa;
    %units_I = stb/day;
    %Pressure of the wells
    P = 275;
    
    P_s = 400;
    I = 1100;
    P0 = 500;
else
    P_b = 400;
    units_Pb = meter^3/day;
    P0 = 100;
end
    % Number of time steps with same pressure (variation in wells)
    tch =2;
    
%% Compute sources
if(use_wells)
    % Set Comp_i=[0,0] in producers to counter X-flow effects...
    
    well_ip = 'ip_tpf';
%     W = verticalWell([], G, rock,  1+ceil(nx/3),   1+ceil(ny/3), [], 'Type', 'bhp', ...
%         'InnerProduct', well_ip, ...
%         'Val', P_s*barsa, 'Radius', 0.125*meter, ...
%         'Name', 'P1', 'Comp_i', [0, 0]);
%     
%     W = verticalWell(W , G, rock,  1+ceil(2*nx/3),   1+ceil(ny/3), [], 'Type', 'bhp', ...
%         'InnerProduct', well_ip, ...
%         'Val', P_s*barsa, 'Radius', 0.125*meter, ...
%         'Name', 'P2', 'Comp_i', [0, 0]);
%     
%     W = verticalWell(W , G, rock,  1+ceil(nx/3),   1+ceil(2*ny/3), [], 'Type', 'bhp', ...
%         'InnerProduct', well_ip, ...
%         'Val', P_s*barsa, 'Radius', 0.125*meter, ...
%         'Name', 'P3', 'Comp_i', [0, 0]);
%     
%     W = verticalWell(W , G, rock, 1+ceil(2*nx/3),   1+ceil(2*ny/3), [], 'Type', 'bhp', ...
%         'InnerProduct', well_ip, ...
%         'Val', P_s*barsa, 'Radius', 0.125*meter, ...
%         'Name', 'P4', 'Comp_i', [0, 0]);
%     
%     W = verticalWell(W , G, rock, ceil(nx/2), ceil(ny/2), [], 'Type', 'bhp',   ...
%         'InnerProduct', well_ip, ...
%         'Val', I*units_I, 'Radius', 0.125*meter, ...
%         'Name', 'I1', 'Comp_i', [1, 0]);
    
    
       W = verticalWell([], G, rock,  1,   1,  [], 'Type', 'bhp', ...
        'InnerProduct', well_ip, ...
        'Val', P_s*barsa, 'Radius', 0.125*meter, ...
        'Name', 'P1', 'Comp_i', [0, 0]);
    
    W = verticalWell(W , G, rock, nx,   1, [], 'Type', 'bhp', ...
        'InnerProduct', well_ip, ...
        'Val', P_s*barsa, 'Radius', 0.125*meter, ...
        'Name', 'P2', 'Comp_i', [0, 0]);
    
    W = verticalWell(W , G, rock,  1,   ny, [], 'Type', 'bhp', ...
        'InnerProduct', well_ip, ...
        'Val', P_s*barsa, 'Radius', 0.125*meter, ...
        'Name', 'P3', 'Comp_i', [0, 0]);
    
    W = verticalWell(W , G, rock, nx,   ny, [], 'Type', 'bhp', ...
        'InnerProduct', well_ip, ...
        'Val', P_s*barsa, 'Radius', 0.125*meter, ...
        'Name', 'P4', 'Comp_i', [0, 0]);
    
    W = verticalWell(W , G, rock, ceil(nx/2), ceil(ny/2), [], 'Type', 'bhp',   ...
        'InnerProduct', well_ip, ...
        'Val', I*units_I, 'Radius', 0.125*meter, ...
        'Name', 'I1', 'Comp_i', [1, 0]);
    
    %% Changing wells parameters
    

    Changing_w
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
x         = initState(G, [], P0*barsa, [0 1]);
x.wellSol = initWellSol(W, 0);
else
%x         = initResSol (G, P0*barsa);
x          = initState(G, [], P0*barsa, [0 1]);
end

%%
linsolve_p = @(S, h) agmg(S, h,  1, 5.0e-11, 1000, 0);
linsolve_t = @(J, F) agmg(J, F, 50, 5.0e-11, 2000, 0);

if(use_wells)
%    psolve = @(x) ...
%       incompTPFA(x, G, T, fluid, 'wells', W, 'LinSolve', linsolve_p);
   tsolve = @(x, dt) ...
      implicitTransport(x, G, dt, rock, fluid, 'wells', W, ...
                        'LinSolve', linsolve_t, 'verbose',false);
else
%    psolve = @(x) ...
%       incompTPFA(x, G, T, fluid, 'bc', bc, 'LinSolve', linsolve_p);
   tsolve = @(x, dt) ...
      implicitTransport(x, G, dt, rock, fluid, 'bc', bc, ...
                        'LinSolve', linsolve_t, 'verbose',false);    
end




%% Create directories to save results
dir = '/home/wagm/cortes/Localdisk/Research/Articles_res/JCP18/02_05/';
%dir = '/dev/media/Sphinx/Doctorado_Delft/articles/2017_CM/Results_1_9/ow/';
%dir = '../2017_CM/Results1/';
%dir = '/dev/media/Sphinx/Doctorado_Delft/articles/2017_CM/Results/';
if(use_wells)
    if(model_SPE)
        folder=[ 'SPE10_' num2str(numel(layers))  'DT_' num2str(DT/day) 'step_' num2str(nstep) 'P_' num2str(P)];
    else
        folder=[ 'per_' num2str(per) 'sz_' num2str(nx)  'layers_' num2str(nz)  'DT_' num2str(DT/day) 'step_' num2str(nstep) 'P_' num2str(P)];
    end
else
    if(model_SPE)
        folder=[ 'SPE10_' num2str(numel(layers))  'DT_' num2str(DT/day) 'step_' num2str(nstep) 'bc_' num2str(P_b)];
    else
        folder=[ 'per_' num2str(per) 'sz_' num2str(nx) '_layers_' num2str(nz)  'DT_' num2str(DT/day) 'step_' num2str(nstep) 'bc_' num2str(P_b)];
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
