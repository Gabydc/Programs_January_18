function[psolve]= ICCG_w(W,tol_p,maxIter,p0,G,hT,fluid,MatOut)
% We initialize the ICCG solver
solver = ICCGSolverAD('tolerance', tol_p,'maxIterations',  maxIter,'x0',p0,'W', W);
fn = @(A, b) solver.solveLinearSystem(A, b);
psolve = @(x) incompTPFA_Def(x, G, hT, fluid, 'wells', W,'LinSolve', fn,'MatrixOutput',MatOut);
end