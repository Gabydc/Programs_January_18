%% Simulation of a Mega-cell Model
% The purpose of this example is to show how one can setup MRST so that it
% is capable of solving models with multi-million cells. To this end, we
% will use the two-point flux-approximation scheme combined with a highly
% efficient algebraic multigrid method (AGMG) to solve the single-phase
% pressure equation
%
% $$\nabla\cdot v = q, \qquad v=\textbf{--}\frac{K}{\mu}\nabla p,$$
%
% for a flow driven by Dirichlet and Neumann boundary conditions. Our
% geological model will be simple a Cartesian grid with anisotropic,
% homogeneous permeability.
close all
clear all
clc
if isempty(mrstPath('agmg'))
   error('This Example Requires the AGMG Module');
end

mrstModule add incomp

%% Define geometry
% Construct a Cartesian grid of size 200-by-100-by-100 (=2,000,000) cells
% in a box of dimensions 10-by-10-by-4 metres. We have been able to run
% this problem with a peak memory use of slightly less than 4GB.
nx = 20; ny = 10; nz = 10;
G = cartGrid([nx, ny, nz], [10, 10, 4]);
s0 = 1;
s  = 1;

%% Process geometry
% The computation of geometric information (centroids and volumes of the
% cells and centroids, normals, and areas for the faces) was accelerated
% significantly in MRST R2016a, but is still faster if one uses the
% C-accelerated routine instead of the standard call to computeGeometry(G)
mrstModule add libgeometry incomp
G = mcomputeGeometry(G);

%% Set rock and fluid data
% The only parameters in the single-phase pressure equation are the
% permeability $K$ and the fluid viscosity $\mu$.
rock = makeRock(G, [1000, 100, 10].* milli*darcy(), 1);
T         = computeTrans(G, rock);

gravity off;
fluid     = initSingleFluid('mu' ,    1*centi*poise     , ...
                            'rho', 1014*kilogram/meter^3);

%% Set boundary conditions
% Here, we impose Neumann conditions (flux of 1 m^3/day) on the global
% left-hand side. Similarly, we set Dirichlet boundary conditions p = 0 on
% the global right-hand side of the grid. For a single-phase
% flow, we need not specify the saturation at inflow boundaries and the
% fluid composition over outflow faces (here, right) is ignored by pside.
resSol = initResSol(G, 0.0);
bc     = fluxside([], G, 'LEFT',  1*meter^3/day());
bc     = pside   (bc, G, 'RIGHT', 0);


%% Solve the linear system
% Solve linear system construced from T and bc to obtain solution for flow
% and pressure in the reservoir. When working with large models, one cannot
% use the standard MLDIVIDE ('\') solver in MATLAB. Here, we use the AGMG
% algebraic multigrid solver.
solver_agmg = false;
if solver_agmg
mrstModule add agmg
resSol = incompTPFA(resSol, G, T, fluid, 'bc', bc, ...
    'LinSolve', @(A,b) agmg(A,b,1));
else
maxIter = 200;
tol =1e-8;
p0 = resSol.pressure;

solver = ICCGSolverAD('tolerance', tol,'maxIterations',  maxIter,'x0',p0,'bc', bc);
linsolve_p = @(A, b) solver.solveLinearSystem(A, b);
psolve = @(x) incompTPFA_Def(x, G, T, fluid, 'bc', bc,'LinSolve', linsolve_p,'MatrixOutput',true);
[resSol,preport]= psolve(resSol);
end
pause
Z = (resSol.A);
solver = DICCGSolverAD('tolerance', tol,'maxIterations',  maxIter,'Z',Z,'x0',p0,'bc', bc);
linsolve_p = @(A, b) solver.solveLinearSystem(A, b);
psolve = @(x) incompTPFA_Def(x, G, T, fluid, 'bc', bc,'LinSolve', linsolve_p);
[resSol,preport]= psolve(resSol);
% clf
% plotCellData(G, convertTo(resSol.pressure(1:G.cells.num), barsa()),'EdgeColor','none');
% title('Cell Pressure [bar]')
% xlabel('x'), ylabel('y'), zlabel('Depth');
% view(3); camproj perspective; axis tight;
% colorbar
%%
 %%
dir ='/dev/media/Sphinx/Doctorado_Delft/2018/February/03_02/';
% % 1. Marker, 2. Color, 3. Markersize, 4. Linestyle, 5. LineWidth
% % linestyles = {'none', '-',  '--', '-.', ':' };
% % markerstyles = {'none', '.', 'o', 'd', '*' '+' 'x' 's'} ;
% % colors = {'r',  [0 0.7 0], 'b', 'm', [0 0.7 0.7], [0.7 0 0.7], [0.5 0.5 0.7]};
% % time = (DT:DT:t)/day;
% optp = { [2 2 3 4 5 6 7], [1 2], 5, [1 2 3 4], [1.8] };
% Plot_val(preport.extras.DA,'savename','Eigenvalues','optp',optp,'x_lab',...
%     'Eigs','y_lab', 'log(eig)', 'revx', true,...
%      'dir', [], 'figure', 10, 'titl', 'Relperm','o_ax', 2);

    %'o_legend', [{'kr_1', 'kr_2'}], ...
    
% Plot_val('dnames',{{'Eigenvalues'}},'optp',optp,'x_lab',...
%     'S','y_lab', 'Permeability', 'dir', dir,'figure', 1);
%%

figure
plotCellData(G, convertTo(resSol.pressure(1:G.cells.num), barsa()), ...
    'EdgeAlpha', 0.050, 'FaceAlpha', 0.375);
title('Cell Pressure [bar]')
xlabel('x'), ylabel('y'), zlabel('Depth');
view(3); camproj perspective; axis tight;
colorbar
%%

Plot_mesh(convertTo(resSol.pressure(1:G.cells.num), barsa()),...
    G, 'dnames', {{'Pressure'}},'o_daspect',[6 6 1], ...
    'l_cb', false, 'o_cb', 1, 'figure', 3, 'showgrid', 2, ...
    'EA', 0.050, 'FA', 0.375, 'titl', 'Cell Pressure [bar]','x_lab', 'x', 'y_lab', 'y', ...
    'inter',   [], 'o_view', [35 35],'o_daspect', [3 3 1.5],  'o_cb', 2,...
    'l_cb', [],  'o_ax',  2, 'dir', [], '3d', true)

%%
T=1;
%(nf,clim,ts,Np,G,nz,time,x)
clear y
clear xy
for i = 1 : 4 
y(:,i) = resSol.pressure/barsa;
xy(i,1) = T;
 %T = (i+1)*ceil(k/Np);
end
Plot_mesh_sub(y, G, 'dnames', {{'Pressure4'}}, ...
    'o_ax', 1, 'l_cb', false, 'o_cb', 1, 'figure', 5, 'showgrid', 2, ...
    'EA', 0, 'FA', 1, 'titl', '','x_lab', '', 'y_lab', '', ...
    'inter',   [], 'o_view', [30 30],'o_daspect', [1 1 1],...
    'l_cb', [],  'dir', [], 'x',xy, 'zd', true)
%%
Plot_mesh_sub(y, G, 'dnames', {{'Pressure4'}}, ...
    'o_ax', 1, 'l_cb', false, 'o_cb', 1, 'figure', 5, 'showgrid', 2, ...
    'EA', 0, 'FA', 1, 'titl', '','x_lab', '', 'y_lab', '', ...
    'inter',   [], 'o_view', [0 90],'o_daspect', [1 1 1],...
    'l_cb', [],  'dir', [], 'x',xy)


%% Copyright notice

% <html>
% <p><font size="-1">
% Copyright 2009-2017 SINTEF ICT, Applied Mathematics.
% </font></p>
% <p><font size="-1">
% This file is part of The MATLAB Reservoir Simulation Toolbox (MRST).
% </font></p>
% <p><font size="-1">
% MRST is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% </font></p>
% <p><font size="-1">
% MRST is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% </font></p>
% <p><font size="-1">
% You should have received a copy of the GNU General Public License
% along with MRST.  If not, see
% <a href="http://www.gnu.org/licenses/">http://www.gnu.org/licenses</a>.
% </font></p>
% </html>
