


if x_true
    %
    np = 8;
else
    np =2;
end

plotv = cell(2,np);

ylab{1} = '||M^{-1}r^k||_2';
ylab{2} = '||M^{-1}r^k||_2/||M^{-1}b||_2';
Name{1} = ['Residual_mr'];
Name{2} = ['Relative_Residual_mr'];
Title{1} = ['Residual ' ylab{1}];
Title{2} = ['Relative Residual ' ylab{2}];

plotv{1,1} = resulte.res;
plotv{2,1} = resulte3.res;
plotv{1,2} = resulte.res / norm(resulte.Mb);
plotv{2,2} = resulte3.res / norm(resulte3.Mb);

if x_true
    
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
    plotv{1,3} = resulte.tres;
    plotv{2,3} = resulte3.tres;
    plotv{1,4} = resulte.tres / norm(resulte.Mb);
    plotv{2,4} = resulte3.tres / norm(resulte3.Mb);
    
    plotv{1,5} = resulte.tresm;
    plotv{2,5} = resulte3.tresm;
    plotv{1,6} = resulte.tresm / norm(resulte.Mb);
    plotv{2,6} = resulte3.tresm / norm(resulte3.Mb);
    
    plotv{1,7} = resulte.terr;
    plotv{2,7} = resulte3.terr;
    plotv{1,8} = resulte.terr / norm(resulte.xtrue);
    plotv{2,8} = resulte3.terr / norm(resulte3.xtrue);  
end






%%
% % 1. Marker, 2. Color, 3. Markersize, 4. Linestyle, 5. LineWidth
% % linestyles = {'none', '-',  '--', '-.', ':' };
% % markerstyles = {'none', '.', 'o', 'd', '*' '+' 'x' 's'} ;
% % colors = {'r',  [0 0.7 0], 'b', 'm', [0 0.7 0.7], [0.7 0 0.7], [0.5 0.5 0.7]};

optp = { [7 3 3 4 5 6 7], [1 2], 4, [1 1 1 1], [1.8] };
ni = 2;
y = cell(ni,1);
for i = 1 : np
    for j = 1 : ni
   y{j} = plotv{j,i};
    end
 Plot_val_s(y,'optp',optp,'x_lab',...
     'Iterations','y_lab',ylab{i},'o_legend', [{['ICCG, min = ' num2str(min(y{1,1}))], ['DICCG, min = ' num2str(min(y{2,1}))] }], ...
    'revx', false, 'dir', [], 'figure', i, 'titl', Title{i},'o_ax', 2) 
end
%%
% res = cell(2,1);
% res{1} = resulte.res;
% res{2} = resulte3.res;
% relres = cell(2,1);
% relres{1} = resulte.res / norm(resulte.Mb);
% relres{2} = resulte3.res;
% 
% if x_true
% tres = cell(2,1);
% tres{1} = resulte.tres;
% tres{2} = resulte3.tres;
% trelres = cell(2,1);
% trelres{1} = resulte.tres / norm(resulte.Mb);
% trelres{2} = resulte3.tres;
% 
% tresm = cell(2,1);
% tresm{1} = resulte.tresm;
% tresm{2} = resulte3.tresm;
% trelresm = cell(2,1);
% trelresm{1} = resulte.tresm / norm(resulte.Mb);
% trelresm{2} = resulte3.tresm;
% 
% terr = cell(2,1);
% terr{1} = resulte.terr;
% trerr{2} = resulte3.terr;
% trelerr = cell(2,1);
% trelerr{1} = resulte.terr/ norm(resulte.xtrue);
% trelerr{2} = resulte3.terr;
% end
distFig('Rows',2,'Columns',np/2)