%... The MatMol Group (2016)
%... Function....... : Automates the simulation for Pe = 1

%... ---------------- %
%... MatMOL Low-order %
%... ---------------- %

    clear all, close all
    Pe = 1; n = 51; zswitch = 0.54; orderD1 = 1; orderD2 = 2;
    TDFR_MatMOL(Pe,n,zswitch,orderD1,orderD2)

%...
    clear all, close all
    Pe = 1; n = 101; zswitch = 0.54; orderD1 = 1; orderD2 = 2;
    TDFR_MatMOL(Pe,n,zswitch,orderD1,orderD2)

%... 
    clear all, close all
    Pe = 1; n = 251; zswitch = 0.54; orderD1 = 1; orderD2 = 2;
    TDFR_MatMOL(Pe,n,zswitch,orderD1,orderD2)


%... ----------------- %
%... MatMOL High-order %
%... ----------------- %

%...
    clear all, close all
    Pe = 1; n = 51; zswitch = 0.54; orderD1 = 4; orderD2 = 4;
    TDFR_MatMOL(Pe,n,zswitch,orderD1,orderD2)

%...
    clear all, close all
    Pe = 1; n = 101; zswitch = 0.54; orderD1 = 4; orderD2 = 4;
    TDFR_MatMOL(Pe,n,zswitch,orderD1,orderD2)

%...
    clear all, close all
    Pe = 1; n = 251; zswitch = 0.54; orderD1 = 4; orderD2 = 4;
    TDFR_MatMOL(Pe,n,zswitch,orderD1,orderD2)
