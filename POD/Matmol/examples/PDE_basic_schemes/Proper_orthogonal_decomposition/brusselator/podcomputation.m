%... The MatMol Group (2016)
%... Computation of the POD basis
    clear all
    clc

%... Load the snapshots (previously computed using the finite element
%... solution in file main_brusselator_fem.m)
    load brus_sol_fem_501

%... Construct the FEM mass matrix for spatial integration
    [MM] = matfem (z, 'neu', 'neu', 'MM');

%... Call the function matpod
    [pods_u] = matpod(u , MM , 'i' , 200);
    [pods_v] = matpod(v , MM , 'i' , 200);

%... save pod_data_brus pods_u pods_v