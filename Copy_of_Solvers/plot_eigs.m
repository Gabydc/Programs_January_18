

np =1;
plotv = cell(2,np);

ylab{1} = 'log(Eigenvalues)';
ylab{2} = 'log(Eigenvectors)';
Name{1} = ['Eigenvalues'];
Name{2} = ['Eigenvectors'];
Title{1} = ['Eigenvalues ' ylab{1}];
Title{2} = ['Eigenvectors ' ylab{2}];

plotv{1,1} = resulte3.DMA;
plotv{2,1} = resulte3.DPMA;

optp = { [7 3 3 4 5 6 7], [1 2], 4, [1 1 1 1], [1.8] };
ni = 2;
y = cell(ni,1);
for i = 1 : np
    for j = 1 : ni
   y{j} = plotv{j,i};
    end
 Plot_val_s(y,'optp',optp,'x_lab',...
     'Iterations','y_lab',ylab{i},'o_legend', [{['ICCG'], ['DICCG'] }], ...
    'revx', true, 'dir', [], 'figure', 1, 'titl', Title{i},'o_ax', 2) 
end
%%
y_lim = { [min(resulte.DMA) max(resulte.DMA)]};
for i = 1 : np
    for j = 1 : ni
   y{j} = plotv{j,i};
    end
 Plot_val_s(y,'optp',optp,'x_lab',...
     'Iterations','y_lab',ylab{i},'o_legend',...
     [{['ICCG, \kappa = ' num2str(resulte.CMA)], ['DICCG, \kappa = '  num2str(resulte3.CPMA)] }], ...
    'revx', true, 'dir', [], 'figure', 2, 'titl', Title{i},'o_ax', 2,'y_lim',y_lim) 
end


%%
%distFig('Rows',2,'Columns',np/2)