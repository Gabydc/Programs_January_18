function varargout = Plot_mesh_sub(y, G, varargin)
%%
% Plot_mesh(convertTo(resSol.pressure(1:G.cells.num), barsa()),...
%     G, 'dnames', {{'Pressure'}},'o_daspect',[6 6 1], ...
%     'l_cb', false, 'o_cb', 1, 'figure', 3, 'showgrid', 2, ...
%     'EA', 0.050, 'FA', 0.375, 'titl', 'Cell Pressure [bar]','x_lab', 'x', 'y_lab', 'y', ...
%     'inter',   [], 'o_view', [35 35],'o_daspect', [3 3 1.5],  'o_cb', 2,...
%     'l_cb', [],  'o_ax',  2, 'dir', [], 'zd',false)
 linestyles = {'none', '-',  '--', '-.', ':' };
%xp = [0.13 0.45 0.13 0.45 0.75];
xp = [0.16 0.48 0.16 0.48 0.8];
xw = [0.6 0.6 0.1 0.1 0.1];      %yp
yp = [0.34 0.34 0.34 0.34 0.05]; %xw
yw = [0.34 0.34 0.34 0.34 0.75];  %yl
%  showgrid = 2;
% % figure
% clf
% plotCellData(G, convertTo(resSol.pressure(1:G.cells.num), barsa()), ...
%      'EdgeAlpha', 0.050, 'FaceAlpha', 0.375,'LineStyle',linestyles{showgrid});
% title('Cell Pressure [bar]')
% xlabel('x'), ylabel('y'), zlabel('Depth');
% view(3); camproj perspective; axis tight;
% colorbar


%%

opt = struct('figure',      1, ...
    'dnames', {{}}, ...
    'showgrid', 1, ...
    'x', [], ...
    'EA', 0, ...
    'FA', 1, ...
    'W',  [], ...
    'titl',   [], ...
    'x_lab',  [], ...
    'y_lab',  [], ...
    'z_lab',  [], ...
    'inter',   [], ...
    'o_view',    [45 45], ...
    'o_daspect', [3 3 3], ....
    'o_cb', [], ...
    'l_cb', [], ...
    'o_ax',  1, ...
    'dir', [], ...
    'zd', false);
[opt] = merge_options(opt, varargin{:});
nf           = opt.figure;
x            = opt.x;
dnames       = opt.dnames;
W            = opt.W;
titl        = opt.titl;
x_lab        = opt.x_lab;
y_lab        = opt.y_lab;
z_lab        = opt.z_lab;
inter        = opt.inter;
o_view       = opt.o_view; ...
o_daspect    = opt.o_daspect;
o_cb         = opt.o_cb;
l_cb         = opt.l_cb;
o_ax         = opt.o_ax;
dir          = opt.dir;
showgrid     = opt.showgrid;
EA           = opt.EA;
FA           = opt.FA;
zd           = opt.zd;


if(zd)
 xp(5) = 0.8;   
end
  clf  
figure(nf);

file = [dnames{1}{1}];
for T = 1 : 4 
ax(T)=axes('position', [xp(T) xw(T) yp(T) yw(T)]);

%ax(T)=subplot(2,2,T);
plotCellData(G, y(:,T),'LineStyle',linestyles{showgrid},'EdgeAlpha', EA, 'FaceAlpha', FA);
%if nz > 1; view(-25,20); end
if (~isempty(x))  
heading = [sprintf(' T = %d days',x(T))]
title([ heading]);
%set(th, 'Position', [xp(T) xw(T) yp(T) yw(T)])
end
if(~zd)
set(get(gca,'title'),'Position',[5.3 10.4 -1.00011])
else
   set(get(gca,'title'),'Position',[5 6.4 -4.000012]) 
end

if o_ax == 1
    axis off
else
    axis equal tight
end
if(~isempty(o_view))
    view(o_view(1),o_view(2))
end
daspect([ o_daspect]);
if (~isempty(W))
    plotWell(G, W(1:4),'color', 'r')
    plotWell(G, W(5), 'color', 'b')
end
if (~isempty(x_lab))
    xlabel(x_lab,'FontSize',16)
end
if (~isempty(y_lab))
    ylabel(y_lab,'FontSize',16)
end
if (~isempty(z_lab))
    ylabel(z_lab,'FontSize',16)
end
if (~isempty(titl))
    title(titl, 'FontSize',16)
end
if (~isempty(inter))
    caxis(inter)
else
    caxis([min(y(:)) max(y(:))])
end
if (o_cb == 2)
    c = colorbar;
else if (o_cb == 3)
        c = colorbar('southoutside');
    else
    end
end
if (l_cb)
    c.Label.String = ['\fontsize{12}' file '[D]'];
end
end
h=colorbar;
set(h, 'Position', [xp(5) xw(5) yp(5) yw(5)])



for i=1:4
      pos=get(ax(i), 'Position');
      set(ax(i), 'Position', [pos(1) 0.87*pos(2) 0.8*pos(3) pos(4)]);
end

if (dir)
    savefigures(gcf, file, dir)
end