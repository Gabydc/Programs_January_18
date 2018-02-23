%... The MatMol Group (2016)
%... Function....... : Automates the simulation for Pe = 10000

%...
    clear all, close all
    Pe = 10000; n = 51; zswitch = 0.54;
    TDFR_MatMOL_StaticRegridding(Pe,n,zswitch)

%...
    clear all, close all
    Pe = 10000; n = 101; zswitch = 0.54;
    TDFR_MatMOL_StaticRegridding(Pe,n,zswitch)

%...
    clear all, close all
    Pe = 10000; n = 251; zswitch = 0.54;
    TDFR_MatMOL_StaticRegridding(Pe,n,zswitch)
