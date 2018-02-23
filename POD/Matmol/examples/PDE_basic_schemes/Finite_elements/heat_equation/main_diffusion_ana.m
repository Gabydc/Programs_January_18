%... The MatMol Group (2016)
%... The diffusion problem solved using separation of variables

    clear all
    clc

%... Spatial coordinates
    nd = 101;                 %... Discretization points
    z  = linspace(0,1,nd)';   

%... Parameters
    k = 0.1;

%... Time span
    tlist = [0:0.05:0.5 0.6:0.1:2 2.5:0.5:4];

%... Analytical solution
    syms xx kk nn

%... Initial conditions
    z0 = 5*(xx^2/2-xx^4/4)+1;
    f1 = z0;
    f2 = z0*(2*kk*cos(nn*pi*xx));
    If1 = int(f1,0,1);
    If2 = int(f2,0,1);

%... Analitical solution
    nbas = 4;
    It = 0;
    for ii = 1:nbas
        ii
        kk_val = cos(ii*pi*z)*exp(-k*ii^2*pi^2*tlist);
        eIf2   = subs(If2,{xx,kk,nn},{z,kk_val,ii});
        It      = It+eIf2;
    end
    eIf1  = subs(If1,{xx},{z});
    x_ana = It + eIf1;
    
%... Plot the solution
    mesh(tlist,z, x_ana)
    xlabel('Time')
    ylabel('z')
    zlabel('x(t,z)')