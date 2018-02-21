%% Set up grid and petrophysical data
% We use a Cartesian grid of size nx-by-ny with homogeneous petrophysical
% data: permeability of 100 mD and porosity of 0.2.
nf = 0;
ceigs = 0;
%Select a linear solver 1.AGMG, 2.GMRES, 3.PCG_IC, 4.DPCG_IC, other Backslash
lsolver=3;
% Select transport solver 1. Explicit, 2. Implicit
tsolver=2;
% We define the size of the reservoir and the inicial points of the grid
% sz is the number of cells, is the same for the two dimensions
sz=35;

 
% If we use wells, we can define here the position of the wells
% w1=10*2^pw;
% w2=2*w1;
nxi=1;
nyi=1;
nx=sz;
ny=sz;
nz=1;
% Compute estimated condition number
cn=0;
%layers = 1:85;
%We define the length of the reservoir (physical dimensions)
Lx=sz;
Ly=sz;
gravity on
% We define the maximum of iterations for the liner solver
maxIterations=500;
% We define the tolerance of the linear solver
% k=8;
% tol = 5*10^(-k);
% We can change the permeability one layer is 10*milli*darcy(), the second
% is 10*milli*darcy()*10^(-per)
x = linspace(0, 1, 11) .';
y = linspace(1, 0, 11) .';

pc_form = 'nonwetting';
muw = 1;
muo = 10;
rhow = 1000;
rhoo = 700;
krw = 2;
kro = 2;
cap_scale = 10;
props = constantProperties([   muw,  muo] .* centi*poise, ...
                           [rhow, rhoo] .* kilogram/meter^3);  
[kr, pc]  = tabulatedSatFunc([x, x.^krw, y.^kro, y.*cap_scale*barsa]);
% Number of deflation vectors
dv=10;
podv{1}=[1:10];
podv{2}=[6:10];
   %% Transport loop vars
    T      = 4800*day();
    %T      = 100*day();
    nstep = 1;
    dT     = T/nstep;
    dT = 20*day;
    dTplot = ceil(T/3);  % plot only every 100th day
    %T = 40*day;
 