% Initial pressure
p0 = x.pressure;

    if(use_wells)
        for i=1:numel(W)
            p0(N+i) = 0;
        end
        if(training)
            W = W1{k};
            %  W(5).val = I_P(k);
        else
            W = W0{k};
            %   W(5).val = I_P(k);
        end
    end
    
    
    %   W(5).val = PI_1;
    if(use_DICCG)
        use_DICCG
        if (window)
            if k < dv+1
                use_ICCG
                if(use_wells)
                    solver = ICCGSolverAD('tolerance', tol,'maxIterations',  maxIter,'x0',p0,'W', W);
                    linsolve_p = @(A, b) solver.solveLinearSystem(A, b);
                    psolve = @(x) incompTPFA_Def(x, G, T, fluid, 'wells', W,'LinSolve', linsolve_p);
                else
                    solver = ICCGSolverAD('tolerance', tol,'maxIterations',  maxIter,'x0',p0,'bc', bc);
                    linsolve_p = @(A, b) solver.solveLinearSystem(A, b);
                    psolve = @(x) incompTPFA_Def(x, G, T, fluid, 'bc',bc,'LinSolve', linsolve_p);
                    
                end
                t0 = tic;
                [x,preport(k)]= psolve(x);
                dt_p(k) = toc(t0);
                POD_V(:,k) = x.pressure;
            else  
                [U,S]=PODbasis(POD_V(:,k-dv:k-1));
                Z = U(:,dpod);
                if(use_wells)
                    for i=1:numel(W)
                        Z(N+i,1)=0;
                    end
                    
                end
                if(use_wells)
                    W = W0{k};
                    %  W(5).val = I_P(k);
                    solver = DICCGSolverAD('tolerance', tol,'maxIterations',  maxIter,'Z',Z,'x0',p0,'W', W,'dir',dir2);
                    linsolve_p = @(A, b) solver.solveLinearSystem(A, b);
                    psolve = @(x) incompTPFA_Def(x, G, T, fluid, 'wells', W,'LinSolve', linsolve_p);
                else
                    solver = DICCGSolverAD('tolerance', tol,'maxIterations',  maxIter,'Z',Z,'x0',p0,'bc', bc,'dir',dir2);
                    linsolve_p = @(A, b) solver.solveLinearSystem(A, b);
                    psolve = @(x) incompTPFA_Def(x, G, T, fluid, 'bc', bc,'LinSolve', linsolve_p);
                end
            end
        else
            if(use_wells)
                W = W0{k};
                for i=1:numel(W)
                    Z(N+i,1)=0;
                end
                %  W(5).val = I_P(k);
                solver = DICCGSolverAD('tolerance', tol,'maxIterations',  maxIter,'Z',Z,'x0',p0,'W', W);
                linsolve_p = @(A, b) solver.solveLinearSystem(A, b);
                psolve = @(x) incompTPFA_Def(x, G, T, fluid, 'wells', W,'LinSolve', linsolve_p);
            else
                solver = DICCGSolverAD('tolerance', tol,'maxIterations',  maxIter,'Z',Z,'x0',p0,'bc', bc);
                linsolve_p = @(A, b) solver.solveLinearSystem(A, b);
                psolve = @(x) incompTPFA_Def(x, G, T, fluid, 'bc', bc,'LinSolve', linsolve_p);
            end
        end
        t0 = tic;
        [x,preport(k)]= psolve(x);
        dt_p(k) = toc(t0);
        
        if(window)
            POD_V(:,k) = x.pressure;
        end
        
        
    else if(use_ICCG)
            use_ICCG
            
            if(use_wells)
                if training
                    W = W1{k};
                    %  W(5).val = I_P(k);
                else
                    W = W0{k};
                    %   W(5).val = I_P(k);
                end
                %   W(5).val = PI_1;
                solver = ICCGSolverAD('tolerance', tol,'maxIterations',  maxIter,'x0',p0,'W', W);
                linsolve_p = @(A, b) solver.solveLinearSystem(A, b);
                psolve = @(x) incompTPFA_Def(x, G, T, fluid, 'wells', W,'LinSolve', linsolve_p);
            else
                solver = ICCGSolverAD('tolerance', tol,'maxIterations',  maxIter,'x0',p0,'bc', bc);
                linsolve_p = @(A, b) solver.solveLinearSystem(A, b);
                psolve = @(x) incompTPFA_Def(x, G, T, fluid, 'bc',bc,'LinSolve', linsolve_p);
                
            end
            t0 = tic;
            [x,preport(k)]= psolve(x);
            dt_p(k) = toc(t0);
        else
            use_agmg
            if(use_wells)
                linsolve_p = @(S, h) agmg(S, h,  1, tol, maxIter, 0);
                psolve = @(x) incompTPFA(x, G, T, fluid, 'wells', W, 'LinSolve', linsolve_p);
            else
                linsolve_p = @(S, h) agmg(S, h,  1, tol, maxIter, 0);
                psolve = @(x) incompTPFA(x, G, T, fluid, 'bc', bc, 'LinSolve', linsolve_p);
            end
            
            t0 = tic;
            x = psolve(x);
            dt_p(k) = toc(t0);
            
            
            
        end
    end
    
    fprintf('[%02d]: Pressure:  %12.5f [s]\n', k, dt_p(k));
    t0 = tic;
    for i=1
        x = tsolve(x, DT/1);
    end
    dt_t(k) = toc(t0);
    fprintf('[%02d]: Transport: %12.5f [s]\n', k, dt_t(k));
    
    t = k * DT;
    if(use_wells)
        [wres{:}]        = prodCurves(W, x, fluid);
        [wellSol{k}]     = x.wellSol;
        
        Prod             = append_wres(Prod, t, wres{:});
        Prod1{k}         = Prod;
        for i=1:5
            pw1(i,k)=x(1).wellSol(i).pressure;
        end
        if training
            I_P(k) = wellSol{k}(5).pressure;
        end
    end
    facePressure(:,k)= x.facePressure;
    Pressure(:,k)    = x.pressure;
    fluxes(:,k)      = x.flux;
    Sat(:,k)         = x.s(:,1);