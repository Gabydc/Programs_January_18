%... The MatMol Group (2016)
%... Function....... : Automates the simulation for Pe = 1

%... ---------------- %
%... MatMOL Low-order %
%... ---------------- %

    clear all, close all
    Pe = 1; n = 51; orderD1 = 1; orderD2 = 2;
    TFBBR_MatMOL(Pe,n,orderD1,orderD2)
    
    clear all, close all
    Pe = 1; n = 101; orderD1 = 1; orderD2 = 2;
    TFBBR_MatMOL(Pe,n,orderD1,orderD2)
    
    clear all, close all
    Pe = 1; n = 251; orderD1 = 1; orderD2 = 2;
    TFBBR_MatMOL(Pe,n,orderD1,orderD2)









%
%... ----------------- %
%... MatMOL High-order %
%... ----------------- %

    clear all, close all
    Pe = 1; n = 51; orderD1 = 4; orderD2 = 4;
    TFBBR_MatMOL(Pe,n,orderD1,orderD2)
    
    clear all, close all
    Pe = 1; n = 101; orderD1 = 4; orderD2 = 4;
    TFBBR_MatMOL(Pe,n,orderD1,orderD2)
    
    clear all, close all
    Pe = 1; n = 251; orderD1 = 4; orderD2 = 4;
    TFBBR_MatMOL(Pe,n,orderD1,orderD2)
