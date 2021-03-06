%% Simulate the SPE10 base case
% The simulation uses a mimetic pressure solver and an implicit transport
% solver.

clear, close all hidden

%VarsSPE10_t
per =1;
Varslayers_t

if(use_wells)
    Prod = struct('t'  , []                  , ...
        'vpt', zeros([0, numel(W)]), ...
        'opr', zeros([0, numel(W)]), ...
        'wpr', zeros([0, numel(W)]), ...
        'wc' , zeros([0, numel(W)]));
    
    append_wres = @(x, t, vpt, opr, wpr, wc) ...
        struct('t'  , [x.t  ; t                  ], ...
        'vpt', [x.vpt; reshape(vpt, 1, [])], ...
        'opr', [x.opr; reshape(opr, 1, [])], ...
        'wpr', [x.wpr; reshape(wpr, 1, [])], ...
        'wc' , [x.wc ; reshape(wc , 1, [])]);
    
    wres = cell([1, 4]);
end
% if (~training)
%     filep=['I_P'];
%     filename=[dir1 filep];
% load(filename)
% end






for k = 1 : nstep,
    p0 = x.pressure;
    if(use_wells)
        for i=1:numel(W)
            p0(N+i) = 0;
        end
    end
    if(use_DICCG)
        use_DICCG
       if k < dv+1
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
           POD_V(:,k) = x.pressure;
       else
           
        

            np = dv;
            dpod = [np-dv+1:np];
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
            solver = DICCGSolverAD('tolerance', tol,'maxIterations',  maxIter,'Z',Z,'x0',p0,'W', W);
            linsolve_p = @(A, b) solver.solveLinearSystem(A, b);
            psolve = @(x) incompTPFA_Def(x, G, T, fluid, 'wells', W,'LinSolve', linsolve_p);
        else
            solver = DICCGSolverAD('tolerance', tol,'maxIterations',  maxIter,'Z',Z,'x0',p0,'bc', bc);
            linsolve_p = @(A, b) solver.solveLinearSystem(A, b);
            psolve = @(x) incompTPFA_Def(x, G, T, fluid, 'bc', bc,'LinSolve', linsolve_p);
        end
        t0 = tic;
        [x,preport(k)]= psolve(x);
        dt_p(k) = toc(t0);
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
    
    
end

%%
%% Save results
if save_res
    if(training)
        for i=1:k
            its(i,1)=preport(1,i).iter;
        end
        ttits_t = sum(its);
        save([dir1  'ttits_t.mat'],'ttits_t')
        if(use_wells)
        save([dir1  'I_P.mat'],'I_P')
        end
    else
        %         for i=1:k
        %             its(i,1)=preport(1,i).iter;
        %         end
        %         ttits = sum(its);
        %         save([dir1  'ttits.mat'],'ttits')
        filetx = ['results.txt'];
        saveits(dir1,filetx,use_ICCG,use_DICCG,use_POD,dpod,k,dv,preport,last)
    end
    filews=['workspace'];
    filename=[dir2 filews];
    save(filename)
    %  clearvars -except Pressure dir1 plot_sol
    if training
        filews=['Pressure'];
        filename=[dir1 filews];
        save(filename,'Pressure')
        
    end
end
%%
if plot_sol
    Plot_t
end