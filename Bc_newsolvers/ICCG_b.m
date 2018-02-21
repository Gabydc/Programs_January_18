function[psolve]= ICCG_b(tol_p,maxIter,p0,G,hT,fluid,bc,MatOut)
% We initialize the ICCG solver
solver = ICCGSolverAD('tolerance', tol_p,'maxIterations', maxIter,'x0',p0);
fn = @(A, b) solver.solveLinearSystem(A, b);
psolve = @(x) incompTPFA_Def(x, G, hT, fluid, 'bc', bc,'LinSolve', fn ,'MatrixOutput',MatOut);
end