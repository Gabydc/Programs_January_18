%... The MatMol Group (2016)
%... Function....... : Automates the simulation for Pe = 10000

%... ------------------------- %
%... Sequencing Method + Euler %
%... ------------------------- %

%...
clear all, close all
Pe = 10000; n = 50; zswitch = 0.54;
TDFR_SeqMeth_Euler(Pe,n,zswitch);

%...
clear all, close all
Pe = 10000; n = 100; zswitch = 0.54;
TDFR_SeqMeth_Euler(Pe,n,zswitch);

%...
clear all, close all
Pe = 10000; n = 250; zswitch = 0.54;
TDFR_SeqMeth_Euler(Pe,n,zswitch);



%... ----------------------------- %
%... Sequencing Method: Transition %
%... ----------------------------- %

%...
clear all, close all
Pe = 10000; n = 50; zswitch = 0.54;
TDFR_SeqMeth_Transition(Pe,n,zswitch);

%...
clear all, close all
Pe = 10000; n = 100; zswitch = 0.54;
TDFR_SeqMeth_Transition(Pe,n,zswitch);

%...
clear all, close all
Pe = 10000; n = 250; zswitch = 0.54;
TDFR_SeqMeth_Transition(Pe,n,zswitch);



%... ----------------------------------- %
%... Sequencing Method: Explicit (ode45) %
%... ----------------------------------- %

%...
clear all, close all
Pe = 10000; n = 50; zswitch = 0.54;
TDFR_SeqMeth_Explicit(Pe,n,zswitch);

%...
clear all, close all
Pe = 10000; n = 100; zswitch = 0.54;
TDFR_SeqMeth_Explicit(Pe,n,zswitch);

%...
clear all, close all
Pe = 10000; n = 250; zswitch = 0.54;
TDFR_SeqMeth_Explicit(Pe,n,zswitch);





%... ------------------------------------ %
%... Sequencing Method: Implicit (ode15s) %
%... ------------------------------------ %

%...
clear all, close all
Pe = 10000; n = 50; zswitch = 0.54;
TDFR_SeqMeth_Implicit(Pe,n,zswitch);

%...
clear all, close all
Pe = 10000; n = 100; zswitch = 0.54;
TDFR_SeqMeth_Implicit(Pe,n,zswitch);

%...
clear all, close all
Pe = 10000; n = 250; zswitch = 0.54;
TDFR_SeqMeth_Implicit(Pe,n,zswitch);



%... ------------------------------------------------ %
%... Sequencing Method: Explicit (ode45) + high-order %
%... ------------------------------------------------ %

%...
clear all, close all
Pe = 10000; n = 50; zswitch = 0.54;
TDFR_SeqMeth_Explicit(Pe,n,zswitch,4);

%...
clear all, close all
Pe = 10000; n = 100; zswitch = 0.54;
TDFR_SeqMeth_Explicit(Pe,n,zswitch,4);

%...
clear all, close all
Pe = 10000; n = 250; zswitch = 0.54;
TDFR_SeqMeth_Explicit(Pe,n,zswitch,4);




%... ---------------------------------------------------------------------------- %
%... Sequencing Method: Explicit (ode45) + convection-diffusion-reaction sequence %
%... ---------------------------------------------------------------------------- %

%...
clear all, close all
Pe = 10000; n = 50; zswitch = 0.54;
TDFR_SeqMeth_Explicit(Pe,n,zswitch,2,'CDR');

%...
clear all, close all
Pe = 10000; n = 100; zswitch = 0.54;
TDFR_SeqMeth_Explicit(Pe,n,zswitch,2,'CDR');

%...
clear all, close all
Pe = 10000; n = 250; zswitch = 0.54;
TDFR_SeqMeth_Explicit(Pe,n,zswitch,2,'CDR');
