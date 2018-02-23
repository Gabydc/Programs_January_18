%... The MatMol Group (2016)
%... Function....... : Automates the simulation for Pe = 10000

%...
    clear all, close all
    Pe = 10000; n = 51;
    TFBBR_MatMOL_StaticRegridding(Pe,n)

%...
    clear all, close all
    Pe = 10000; n = 101;
    TFBBR_MatMOL_StaticRegridding(Pe,n)

%...
    clear all, close all
    Pe = 10000; n = 251;
    TFBBR_MatMOL_StaticRegridding(Pe,n)


