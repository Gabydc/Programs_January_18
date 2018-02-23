%... The MatMol Group (2016)
%... Function....... : Automates the simulation for Pe = 100

%... ------------------------- %
%... Sequencing Method + Euler %
%... ------------------------- %

%... Pe = 100
clear all, close all
Pe = 100; n = 50;
TFBBR_SeqMeth_Euler(Pe,n);

%...
clear all, close all
Pe = 100; n = 100;
TFBBR_SeqMeth_Euler(Pe,n);

%...
clear all, close all
Pe = 100; n = 250;
TFBBR_SeqMeth_Euler(Pe,n);





%... ----------------------------- %
%... Sequencing Method: Transition %
%... ----------------------------- %

%... Pe = 100
clear all, close all
Pe = 100; n = 50;
TFBBR_SeqMeth_Transition(Pe,n);

%...
clear all, close all
Pe = 100; n = 100;
TFBBR_SeqMeth_Transition(Pe,n);

%...
clear all, close all
Pe = 100; n = 250;
TFBBR_SeqMeth_Transition(Pe,n);





%... ----------------------------------- %
%... Sequencing Method: Explicit (ode45) %
%... ----------------------------------- %

%...
clear all, close all
Pe = 100; n = 50;
TFBBR_SeqMeth_Explicit(Pe,n);

%...
clear all, close all
Pe = 100; n = 100;
TFBBR_SeqMeth_Explicit(Pe,n);

%...
clear all, close all
Pe = 100; n = 250;
TFBBR_SeqMeth_Explicit(Pe,n);







%... ------------------------------------ %
%... Sequencing Method: Implicit (ode15s) %
%... ------------------------------------ %

%...
clear all, close all
Pe = 100; n = 50;
TFBBR_SeqMeth_Implicit(Pe,n);

%...
clear all, close all
Pe = 100; n = 100;
TFBBR_SeqMeth_Implicit(Pe,n);

%...
clear all, close all
Pe = 100; n = 250;
TFBBR_SeqMeth_Implicit(Pe,n);








%... ------------------------------------------------ %
%... Sequencing Method: Explicit (ode45) + high-order %
%... ------------------------------------------------ %

%...
clear all, close all
Pe = 100; n = 50;
TFBBR_SeqMeth_Explicit(Pe,n,4);

%...
clear all, close all
Pe = 100; n = 100;
TFBBR_SeqMeth_Explicit(Pe,n,4);

%...
clear all, close all
Pe = 100; n = 250;
TFBBR_SeqMeth_Explicit(Pe,n,4);





%... ---------------------------------------------------------------------------- %
%... Sequencing Method: Explicit (ode45) + convection-diffusion-reaction sequence %
%... ---------------------------------------------------------------------------- %

%...
clear all, close all
Pe = 100; n = 50;
TFBBR_SeqMeth_Explicit(Pe,n,2,'CDR');

%...
clear all, close all
Pe = 100; n = 100;
TFBBR_SeqMeth_Explicit(Pe,n,2,'CDR');

%...
clear all, close all
Pe = 100; n = 250;
TFBBR_SeqMeth_Explicit(Pe,n,2,'CDR');
