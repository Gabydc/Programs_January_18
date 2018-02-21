% We initialize the solver
solver = ICCGSolverAD('tolerance', tol_p,'maxIterations', maxIter,'x0',p0);
fn = @(A, b) solver.solveLinearSystem(A, b);
psolve = @(x) incompTPFA_Def(x, G, hT, fluid, 'bc', bc,'LinSolve', fn ,'MatrixOutput',true);