%... The MatMol Group (2016)
%... Function....... : Automates the simulation for Pe = 10000


%... ---------------- %
%... MatMOL Low-order %
%... ---------------- %

    clear all, close all
    Pe = 10000; n = 51; zswitch = 0.54; orderD1 = 1; orderD2 = 2;
    TDFR_MatMOL(Pe,n,zswitch,orderD1,orderD2)
    
    clear all, close all
    Pe = 10000; n = 101; zswitch = 0.54; orderD1 = 1; orderD2 = 2;
    TDFR_MatMOL(Pe,n,zswitch,orderD1,orderD2)
    
    clear all, close all
    Pe = 10000; n = 251; zswitch = 0.54; orderD1 = 1; orderD2 = 2;
    TDFR_MatMOL(Pe,n,zswitch,orderD1,orderD2)
    
    
    
    
     
    
%... ----------------------------------- %
%... MatMOL Low-order + Forcing function %
%... ----------------------------------- %
    
    clear all, close all
    Pe = 10000; n = 51; zswitch = 0.54; orderD1 = 1; orderD2 = 2;
    TDFR_MatMOL_ForcingFunction(Pe,n,zswitch,orderD1,orderD2)
    
    clear all, close all
    Pe = 10000; n = 101; zswitch = 0.54; orderD1 = 1; orderD2 = 2;
    TDFR_MatMOL_ForcingFunction(Pe,n,zswitch,orderD1,orderD2)
    
    clear all, close all
    Pe = 10000; n = 251; zswitch = 0.54; orderD1 = 1; orderD2 = 2;
    TDFR_MatMOL_ForcingFunction(Pe,n,zswitch,orderD1,orderD2)
    
    
    
    
    
    
%... ------------------------------ %
%... MatMOL Low-order + Elimination %
%... ------------------------------ %
    
    clear all, close all
    Pe = 10000; n = 51; zswitch = 0.54; orderD1 = 1; orderD2 = 2;
    TDFR_MatMOL_Elimination(Pe,n,zswitch,orderD1,orderD2)
    
    clear all, close all
    Pe = 10000; n = 101; zswitch = 0.54; orderD1 = 1; orderD2 = 2;
    TDFR_MatMOL_Elimination(Pe,n,zswitch,orderD1,orderD2)
    
    
    clear all, close all
    Pe = 10000; n = 251; zswitch = 0.54; orderD1 = 1; orderD2 = 2;
    TDFR_MatMOL_Elimination(Pe,n,zswitch,orderD1,orderD2)
    
        
    
    
%... ----------------- %
%... MatMOL High-order %
%... ----------------- %
    
    clear all, close all
    Pe = 10000; n = 51; zswitch = 0.54; orderD1 = 4; orderD2 = 4;
    TDFR_MatMOL(Pe,n,zswitch,orderD1,orderD2)
    
    clear all, close all
    Pe = 10000; n = 101; zswitch = 0.54; orderD1 = 4; orderD2 = 4;
    TDFR_MatMOL(Pe,n,zswitch,orderD1,orderD2)
    
    clear all, close all
    Pe = 10000; n = 251; zswitch = 0.54; orderD1 = 4; orderD2 = 4;
    TDFR_MatMOL(Pe,n,zswitch,orderD1,orderD2)
    
    
    

%... ------------------------------------ %
%... MatMOL High-order + Forcing function %
%... ------------------------------------ %
    
    clear all, close all
    Pe = 10000; n = 51; zswitch = 0.54; orderD1 = 4; orderD2 = 4;
    TDFR_MatMOL_ForcingFunction(Pe,n,zswitch,orderD1,orderD2)
    
    clear all, close all
    Pe = 10000; n = 101; zswitch = 0.54; orderD1 = 4; orderD2 = 4;
    TDFR_MatMOL_ForcingFunction(Pe,n,zswitch,orderD1,orderD2)
    
    clear all, close all
    Pe = 10000; n = 251; zswitch = 0.54; orderD1 = 4; orderD2 = 4;
    TDFR_MatMOL_ForcingFunction(Pe,n,zswitch,orderD1,orderD2)
    
    
    
    
%... ------------------------------- %
%... MatMOL High-order + Elimination %
%... ------------------------------- %
    
    clear all, close all
    Pe = 10000; n = 51; zswitch = 0.54; orderD1 = 4; orderD2 = 4;
    TDFR_MatMOL_Elimination(Pe,n,zswitch,orderD1,orderD2)
    
    clear all, close all
    Pe = 10000; n = 101; zswitch = 0.54; orderD1 = 4; orderD2 = 4;
    TDFR_MatMOL_Elimination(Pe,n,zswitch,orderD1,orderD2)
    
    
    clear all, close all
    Pe = 10000; n = 251; zswitch = 0.54; orderD1 = 4; orderD2 = 4;
    TDFR_MatMOL_Elimination(Pe,n,zswitch,orderD1,orderD2)
    
    
    
