%... The MatMol Group (2016)
%... Function....... : Automates the simulation for Pe = 10000


%... --------------------- %
%... MatMOL: Flux Limiter  %
%... --------------------- %

%... Minmod
    clear all, close all
    Pe = 10000; n = 51; orderD1 = 4; orderD2 = 4;
    TFBBR_MatMOL_FluxLimiter(Pe,n,orderD1,orderD2)

%...
    clear all, close all
    Pe = 10000; n = 101; orderD1 = 4; orderD2 = 4;
    TFBBR_MatMOL_FluxLimiter(Pe,n,orderD1,orderD2)

%...
    clear all, close all
    Pe = 10000; n = 251; orderD1 = 4; orderD2 = 4;
    TFBBR_MatMOL_FluxLimiter(Pe,n,orderD1,orderD2)

%... Koren
    clear all, close all
    Pe = 10000; n = 51; orderD1 = 4; orderD2 = 4;
    TFBBR_MatMOL_FluxLimiter(Pe,n,orderD1,orderD2,1.0e-4,1.0e-6,'Koren')

%...
    clear all, close all
    Pe = 10000; n = 101; orderD1 = 4; orderD2 = 4;
    TFBBR_MatMOL_FluxLimiter(Pe,n,orderD1,orderD2,1.0e-4,1.0e-6,'Koren')

%...
    clear all, close all
    Pe = 10000; n = 251; orderD1 = 4; orderD2 = 4;
    TFBBR_MatMOL_FluxLimiter(Pe,n,orderD1,orderD2,1.0e-4,1.0e-6,'Koren')

%... MC
    clear all, close all
    Pe = 10000; n = 51; orderD1 = 4; orderD2 = 4;
    TFBBR_MatMOL_FluxLimiter(Pe,n,orderD1,orderD2,1.0e-4,1.0e-6,'MC')

%...
    clear all, close all
    Pe = 10000; n = 101; orderD1 = 4; orderD2 = 4;
    TFBBR_MatMOL_FluxLimiter(Pe,n,orderD1,orderD2,1.0e-4,1.0e-6,'MC')

%...
    clear all, close all
    Pe = 10000; n = 251; orderD1 = 4; orderD2 = 4;
    TFBBR_MatMOL_FluxLimiter(Pe,n,orderD1,orderD2,1.0e-4,1.0e-6,'MC')

%... Smart
    clear all, close all
    Pe = 10000; n = 51; orderD1 = 4; orderD2 = 4;
    TFBBR_MatMOL_FluxLimiter(Pe,n,orderD1,orderD2,1.0e-4,1.0e-6,'Smart')

%...
    clear all, close all
    Pe = 10000; n = 101; orderD1 = 4; orderD2 = 4;
    TFBBR_MatMOL_FluxLimiter(Pe,n,orderD1,orderD2,1.0e-4,1.0e-6,'Smart')

%...
    clear all, close all
    Pe = 10000; n = 251; orderD1 = 4; orderD2 = 4;
    TFBBR_MatMOL_FluxLimiter(Pe,n,orderD1,orderD2,1.0e-4,1.0e-6,'Smart')

%... Superbee
    clear all, close all
    Pe = 10000; n = 51; orderD1 = 4; orderD2 = 4;
    TFBBR_MatMOL_FluxLimiter(Pe,n,orderD1,orderD2,1.0e-4,1.0e-6,'Superbee')

%...
    clear all, close all
    Pe = 10000; n = 101; orderD1 = 4; orderD2 = 4;
    TFBBR_MatMOL_FluxLimiter(Pe,n,orderD1,orderD2,1.0e-4,1.0e-6,'Superbee')

%...
    clear all, close all
    Pe = 10000; n = 251; orderD1 = 4; orderD2 = 4;
    TFBBR_MatMOL_FluxLimiter(Pe,n,orderD1,orderD2,1.0e-4,1.0e-6,'Superbee')


%... Van Leer
    clear all, close all
    Pe = 10000; n = 51; orderD1 = 4; orderD2 = 4;
    TFBBR_MatMOL_FluxLimiter(Pe,n,orderD1,orderD2,1.0e-4,1.0e-6,'VanLeer')

%...
    clear all, close all
    Pe = 10000; n = 101; orderD1 = 4; orderD2 = 4;
    TFBBR_MatMOL_FluxLimiter(Pe,n,orderD1,orderD2,1.0e-4,1.0e-6,'VanLeer')

%...
    clear all, close all
    Pe = 10000; n = 251; orderD1 = 4; orderD2 = 4;
    TFBBR_MatMOL_FluxLimiter(Pe,n,orderD1,orderD2,1.0e-4,1.0e-6,'VanLeer')
