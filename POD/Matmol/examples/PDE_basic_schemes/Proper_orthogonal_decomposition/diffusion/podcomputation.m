%... The MatMol Group (2016)
%... Computation of the POD basis
    clear all
    clc

%... Load the snapshots (in order to save time, the snapshots were
%... previously computed using file main_diffusion_ana.m)
    load data_snap_dif

%... Construct the FEM mass matrix for spatial integration using functio
%... matfem.m
    [MM] = matfem (z, 'neu', 'neu', 'MM');

%... Call the function matpod which computes the POD basis
    [pods] = matpod(x_ana , MM , 'i' , 200);