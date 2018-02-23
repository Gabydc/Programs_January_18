%... Computation of the POD basis from the data of the spiral
%... and the front 
    clear all
    clc

%... Load the FEM matrices (previously computed for this problem to save time)
    load data/matrix_fem_1925
    nd = size(MM,1);
    clear fem DM iMM fem ms

%... Load the snapshots (previously computed for this problem to save time using file main_fhn.m)
    load data/front_data_201311.mat   % Front snapshots
    vf = vv;
    wf = ww;
    clear vv ww
    load data/spiral_data_201311.mat  % Spiral snapshots
    vv = [vf vv];
    ww = [wf ww];
    clear vf wf

%... Compute the PODs
    pods_v   = matpod(vv , MM , 'i' , 100);
    pods_w   = matpod(ww , MM , 'i' , 100);
    phi_v    = pods_v.phi;
    phi_w    = pods_w.phi;
    lambda_v = diag(pods_v.lambda);
    lambda_w = diag(pods_w.lambda);