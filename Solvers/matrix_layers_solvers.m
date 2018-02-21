clear all
close all
clc
x = 3;
y = 3;
l = 3;
s0 = 1;
s  = 1;
nel  = x*y*l-2*x;
[A,b,Z] = matrixf(x,y,l,s0,s);


solver.x0(size(b,1),1) = rand;
solver.maxIterations = 1500;
solver.tolerance = 5*10^-5;

%% Solve the matrix with DICCG and DICCG

[L] = ichol(A);
%
%  % Display message of options
%  dopts      = opts{1};
%  % Compute condition number of the matrix A
%  A_cn       = opts{2};
%  % Compute true solution
%  x_true     = opts{3};
%  % Checks the convergence of the method, in the case of DICCG the residual
%  % can increas, in that case, the solution will be the solution with
%  % minimal residual
%  Convergence = opts{4};
%  % Save the variables to plot residual
%  Residual    = opts{5};
%  % Display the number of iterations
%  Iter_m      = opts{6};
%  % Compute eigenvalues of matrix A
%  Amatrix_eigs = opts{7};
%  % Compute eigenvalues of matrix M^{-1}A
%  MAmatrix_eigs = opts{8};
 x_true  = true;
disp('ICCG')
opts = {{true, false, x_true , true, true, true, true, true}};
[result,flag,res,its,resvec,resulte]  = ICCG_MRST(A,b,...
    solver.tolerance,min(solver.maxIterations,nel),L,L',solver.x0, ...
    'opts',opts);



opts = {{true, false, x_true, true, true, true, true, true, true,true}};
disp('DICCG')
[result3,flag3,res3,its3,resvec3,resulte3]  = DICCG_MRST(A,b,Z,...
    solver.tolerance,nel,L,L',solver.x0, ...
    'opts',opts);
%%
% plot_results
%plot_eigs
 
