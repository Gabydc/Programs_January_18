%% Simulate the SPE10 base case
% The simulation uses a mimetic pressure solver and an implicit transport
% solver.

clear, close all hidden

%VarsSPE10_t
VarsSPE_sp

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
% if (use_DICCG) && (~window)
%     filep=['Pressure'];
%     filename=[dir1 filep];
%     load(filename)
%     [U,S]=PODbasis(Pressure);
%     Z = U(:,dpod);
% end






for k = 1 : numel(W)-1
    
    p0 = x.pressure;        
       W = W1{k};
        %W.val
    if(use_wells)
        for i=1:numel(W)
            p0(N+i) = 0;
        end

 
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
        
        fprintf('[%02d]: Pressure:  %12.5f [s]\n', k, dt_p(k));
        
        t = k ;
         if(use_wells)
             [wres{:}]        = prodCurves(W, x, fluid);
             [wellSol{k}]     = x.wellSol;
%             
%             Prod             = append_wres(Prod, t, wres{:});
%             Prod1{k}         = Prod;
            for i=1:numel(W)
                pw1(i,k)=x(1).wellSol(i).pressure;
            end
%             if training
%                 I_P(k) = wellSol{k}(5).pressure;
%             end
         end
         facePressure(:,k)= x.facePressure;
         Pressure(:,k)    = x.pressure;
         fluxes(:,k)      = x.flux;
         Sat(:,k)         = x.s(:,1);
        figure
     plotCellData(G, POD_V(:,k)/barsa,'LineStyle','none');
   
        
end
%%
 W = W0{1};
k = k + 1;
        t0 = tic;
        [x,preport(k)]= psolve(x);
        dt_p(k) = toc(t0);
        POD_V(:,k) = x.pressure;
        
[U,S]=PODbasis(POD_V(:,1:5));
%Z = U(:,dpod);
Z=POD_V(:,1:4);
nf=0;
if use_DICCG
    nf = nf+1;
    file{nf} = ['eig_pod'];
    f(nf) = figure(nf);
    plot(log((diag(S))),'*r');
    set(gca, 'XDir','reverse')
    ylabel('log(Value) ','FontSize',16)
    xlabel('Eigenvalue','FontSize',16)
    axis('tight');
    nf = nf + 1;
    file{nf} = ['eig_vect' ];
    f(nf) = figure(nf);
    plot(Z);
    ylabel('vectors ','FontSize',16)
    xlabel('eigenvectors','FontSize',16)
    axis('tight')
end
pause
    for i=1:numel(W)
        Z(N+i,1)=0; 
        p0(N+i) = 0;
    end

if(use_wells)
    W = W0{1};
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

facePressure(:,k)= x.facePressure;
Pressure(:,k)    = x.pressure;
fluxes(:,k)      = x.flux;
Sat(:,k)         = x.s(:,1);
[wres{:}]        = prodCurves(W, x, fluid);
[wellSol{k}]     = x.wellSol;

for i=1:numel(W)
    pw1(i,k)=x(1).wellSol(k).pressure;
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
        filetx = ['results.txt'];
        
        if(use_DICCG)
            saveits_w(dir1,filetx,use_ICCG,use_DICCG,use_POD,dpod,k,dv,preport,last)
        else
            saveits(dir1,filetx,use_ICCG,use_DICCG,use_POD,dpod,k,dv,preport,last)
        end
    else
        %         for i=1:k
        %             its(i,1)=preport(1,i).iter;
        %         end
        %         ttits = sum(its);
        %         save([dir1  'ttits.mat'],'ttits')
        filetx = ['results.txt'];
        % Save only DICCG iterations
        %saveits(dir1,filetx,use_ICCG,use_DICCG,use_POD,dpod,k,dv,preport,last)
        %Save ICCG and DICCG
        saveits_w(dir1,filetx,use_ICCG,use_DICCG,use_POD,dpod,k,dv,preport,last)
        
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
DT = 1;
nstep = k;
t = k;
if plot_sol
    Plot_t1
end