function varargout = Plot_1(x, G, varargin)
opt = struct('field',       'bhp', ...
    'timescale',   'days', ...
    'figure',      0, ...
    'dnames', {{}}, ...
    'DT', [], ...
    't',  [], ...
    'W',  [], ...
    'title',   [], ...
    'xlabel',  [], ...
    'ylabel',  [], ...
    'o_view',    [45 45], ...
    'o_daspect', [3 3 3], ....
    'o_cb', [], ...
    'l_cb', [], ...
    'o_ax',  1, ...
    'o_save', []);
[opt] = merge_options(opt, varargin{:});
lowermargin  = opt.lowermargin;
plotwidth    = opt.plotwidth;
linewidth    = opt.linewidth;
field        = opt.field;
linestyles   = opt.linestyles;
markerstyles = opt.markerstyles;
pcol         = opt.colors;
timescale    = opt.timescale;
nf           = opt.figure;
dnames = opt.dnames;
DT           = opt.DT;
t            = opt.t;
W            = opt.W;
title        = opt.title;
xlabel       = opt.xlabel;
ylabel       = opt.ylabel;
o_view       = opt.o_view; ...
o_daspect    = opt.o_daspect;
o_cb         = opt.o_cb;
l_cb         = opt.l_cb;
o_ax         = opt.o_ax;
o_save       = opt.o_save;
clf

if (DT)
    time = (DT:DT:t)/day;
else
    time = 1;
end
figure(nf);
file = [dnames{1}{1}];
plotCellData(G, x,'LineStyle',linestyles{1});
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
caxis([min(x) max(x)])
if (o_cb == 1)
    c = colorbar;
else if (o_cb == 2)
        c = colorbar('southoutside');
    else
    end
end
if (l_cb)
    c.Label.String = ['\fontsize{12}' file '[D]'];
end
if (o_save)
        savefigures(gcf, file, o_save)
end