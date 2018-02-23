%%
% Two-phase Incompressible Flow Solver - Gabriela Diaz 17/05/17
% This program simulates two-phase flow through a porous media. The
% pressure equation is solved using DICCG Method and is compared with ICCG
close all
clear all
clc
% The capillary pressure can be taked into account (cp =1)
for cp = [0 ]
    
    % The contrast between permeability layers can be varied one layer has
    % permeability of 10*milli*darcy, the secon is varied as
    % 10*10^(-per)*milli*darcy
    for per=[1]
        
        % In this part the solver is chosen 0 means ICCG (no deflation) and
        % 1 is with deflation
        for def=[0 1]
            
           for tol = [ 8 ]
                    % Initialize the variables, size, fluid values, time
                    % step, ...
                    varsbc
                    %% Transport loop vars
                    %T      = 4800*day();
                    T      = 133*day();
                    nstep = 11;
                    dT     = T/nstep;
                    %dT = 20*day;
                    dTplot = ceil(T/3);  % plot only every 100th day
                    %T = 40*day;
                    
                    
                    tol_p = 5*10^(-tol);
                    tol_t   = 5.0e-7;
                    maxIter = 1000;
                    dpod =[];
                    pod =0;
                    %Create the directory
                    dir='/mnt/sda2/cortes/Results/2018/02/21/plots/'; 
                    folder=[ '10-' num2str(tol) '_' num2str(sz) 'nz' num2str(nz) 'perm_' num2str(per) 'cp' num2str(cp)];
                    %mkdir([dir], folder)
                    dir1 = [dir folder '/'];       
                    folder=[  'def_' num2str(def) '_pod_' num2str(numel(dpod)) ];
                    %mkdir([dir1], folder)
                    dir2 = [dir1 folder '/'];
                    name = [dir2 'workspace'];
                    
                    load (name)
                    def
%                     pause
                    if def == 0
                        nm =1;
                    else  
  
                            nm = 2;
                    end
                
                        terr{nm,1}  = preport(11).extras.terr;
                        nxtrue(nm)  = norm(preport(11).extras.xtrue);
                        nterr{nm,1} = terr{nm,1} / nxtrue(nm); 
                        tres{nm,1}  = preport(11).extras.tres;
                        nb(nm)    = norm(preport(11).extras.b);
                        ntres{nm,1}  = tres{nm,1} / nb(nm);
                        tresm{nm,1} = preport(11).extras.tresm;
                        nresm(nm)   = norm(preport(11).extras.b);
                        ntresm{nm,1}  = tresm{nm,1} / nresm(nm);
                        res{nm,1} = preport(11).extras.res;
                        nmb(nm)   = norm(preport(11).extras.Mb);
                        nres{nm,1}  = res{nm,1} / nmb(nm);
                       
                    clearvars -except per def pod cp terr tres tresm  nterr ntres ntresm nres tol 
                end
        end
    end
end
%%

ylab{1} = '||M^{-1}r^k||_2';
ylab{2} = '||M^{-1}r^k||_2/||M^{-1}b||_2';
Name{1} = ['Residual_mr'];
Name{2} = ['Relative_Residual_mr'];
Title{1} = ['Residual ' ylab{1}];
Title{2} = ['Relative Residual ' ylab{2}];

    ylab{3} = '||M^{-1}(b-Ax^k)||_2';
    ylab{4} = '||M^{-1}(b-Ax^k)||_2/||M^{-1}b||_2';
    ylab{5} = '||b-Ax^k||_2';
    ylab{6} = '||b-Ax^k||_2/||b||_2';
    ylab{7} = '||x-x^k||_2';
    ylab{8} = '||x-x^k||_2/||x||_2';
    
    Name{3} = ['True_Residual_r'];
    Name{4} = ['True_Relative_Residual_mr'];
    Name{5} = ['True_Residual_r'];
    Name{6} = ['True_Relative_Residual_r'];
    Name{7} = ['True_Error'];
    Name{8} = ['True_Relative_Error'];
    
    Title{3} = ['True Residual ' ylab{3}];
    Title{4} = ['True Relative Residual ' ylab{4}];
    Title{5} = ['True Residual ' ylab{5}];
    Title{6} = ['True Relative Residual ' ylab{6} ];
    Title{7} = ['True Error ' ylab{7}];
    Title{8} = ['True Relative Error ' ylab{8}];

%linestyles = {'none', '-',  '--', '-.', ':' };
%markerstyles = {'none', '.', 'o', 'd', '*' '+' 'x' 's'} ;
%colors = {'r',  [0 0.7 0], 'b', 'm', [0 0.7 0.7], [0.7 0 0.7], [0.5 0.5 0.7]};
% optp is a vector containing the options to plot
% the elements are arranged as follows:
% 1. Marker, 2. Color, 3. Markersize, 4. Linestyle, 5. LineWidth
% Example
optp = { [3 7 6 4 5 6 7], [5 4], 5, [1 1 1 1], [1.8] };
val{2} = nres;
val{4} = ntresm; 
val{6} = ntres; 
val{8} = nterr; 
for i = [ 2 4 6 8]
    folder=[ '/mnt/sda2/cortes/Results/2018/02/21/plots/6_tol' num2str(tol) '/p' num2str(per) '/' ];
                    mkdir( folder)
                    dir11 = [ folder ];
 Plot_val_s(val{i},'optp',optp,'x_lab',...
     'Iterations','y_lab',ylab{i},'o_legend', [{['ICCG'], ['DICCG' ] }], ...
    'revx', false, 'dir', dir11, ...
    'figure', i, 'titl', Title{i},'o_ax', 2,'savename',Name{i}) 
end
