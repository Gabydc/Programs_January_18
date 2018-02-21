function varargout = Plot_2(y, varargin)
%%
linestyles = {'none', '-',  '--', '-.', ':' };
markerstyles = {'.', 'o', 'd', '*' '+' 'x' 's'} ;
colors = {'r', 'm' [0 0.7 0], 'b',  [0 0.7 0.7], [0.7 0 0.7], [0.5 0.5 0.7]};
% optp is a vector containing the options to plot
% the elements are arranged as follows:
% 1. Marker, 2. Color, 3. Markersize, 4. Linestyle, 5. LineWidth
%% Example
% optp = [ 1 2 8 1 1 ];
% x = 1:100;
% y = sin(2*x);
% figure
% plot(x,y,'Marker',markerstyles{optp(1)},'MarkerFaceColor', ...
%     colors{optp(2)}, 'MarkerSize',optp(3), 'Linestyle', ...
%     linestyles{optp(4)},'linewidth', optp(5));

    %%
opt = struct('optp',   [ 1 2 8 1 1 ], ...
    'field',       'bhp', ...
    'timescale',   'days', ...
    'figure',      0, ...
    'dnames', {{}}, ...
    'DT', [], ...
    't',  [], ...
    'W',  [], ...
    'titl',   [], ...
    'x_lab',  [], ...
    'y_lab',  [], ...
    'o_view',    1, ...
    'o_daspect', [3 3 3], ....
    'o_cb', [], ...
    'l_cb', [], ...
    'o_ax',  2, ...
    'x', [], ...
    'o_legend','', ...
    'dir', []);
[opt] = merge_options(opt, varargin{:});
optp         = opt.optp;
field        = opt.field;
timescale    = opt.timescale;
nf           = opt.figure;
dnames = opt.dnames;
DT           = opt.DT;
t            = opt.t;
W            = opt.W;
titl         = opt.titl;
x_lab        = opt.x_lab;
y_lab        = opt.y_lab;
o_view       = opt.o_view; ...
o_daspect    = opt.o_daspect;
o_cb         = opt.o_cb;
l_cb         = opt.l_cb;
o_ax         = opt.o_ax;
x            = opt.x;
o_legend     = opt.o_legend;
dir          = opt.dir;

if (DT)
    time = (DT:DT:t)/day;
else
    time = 1;
end
nf = nf + 1;
figure(nf);
file = [dnames{1}{1}];
if (x)    
f = plot(x, y,'Marker',markerstyles{optp(1)},'MarkerFaceColor', ...
    colors{optp(2)}, 'MarkerSize',optp(3), 'Linestyle', ...
    linestyles{optp(4)},'linewidth', optp(5));
else     
f = plot(y,'Marker',markerstyles{optp(1)},'MarkerFaceColor', ...
    colors{optp(2)}, 'MarkerSize',optp(3), 'Linestyle', ...
    linestyles{optp(4)},'linewidth', optp(5));   
end

if (isempty(o_legend))
else
legend(o_legend)
end
if (~isempty(x_lab))
xlabel(x_lab,'FontSize',16)
end
if (~isempty(y_lab))
ylabel(y_lab,'FontSize',16)
end
if (~isempty(titl))
title(titl, 'FontSize',16)
end
    
if o_ax == 1
    axis off
else if o_ax == 2
    axis tight
    else
     axis off equal tight
    end
end

if (o_cb == 1)
    c = colorbar;
else if (o_cb == 2)
        c = colorbar('southoutside');
    end
end
if (l_cb)
    c.Label.String = ['\fontsize{12}' file '[D]'];
end
if (~isempty(dir))
        savefigures(gcf, file, dir)
end