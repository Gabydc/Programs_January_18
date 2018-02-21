function[fn]= ICCG_b(tol_p,maxIter,p0)
% We initialize the ICCG solver
solver = ICCGSolverAD('tolerance', tol_p,'maxIterations', maxIter,'x0',p0);
fn = @(A, b) solver.solveLinearSystem(A, b);
end