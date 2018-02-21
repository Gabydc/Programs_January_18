



    switch lsolver
        case 1
            mrstModule add agmg
            solver = AGMGSolverAD('tolerance', tol_p);
        case 2
            solver = GMRES_ILUSolverAD('tolerance', tol_p);
        case 3
            solver = ICCGSolverAD('tolerance', tol_p,'maxIterations', maxIter,'x0',p0);
        case 4
            solver = DICCGSolverAD('tolerance', tol_p,'maxIterations', maxIter,'Z',Z,'x0',p0);
        otherwise
            solver = BackslashSolverAD();
    end