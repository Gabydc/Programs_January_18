function varargout = Plot_val(y, varargin)
%%
% % 1. Marker, 2. Color, 3. Markersize, 4. Linestyle, 5. LineWidth
% % linestyles = {'none', '-',  '--', '-.', ':' };
% % markerstyles = {'none', '.', 'o', 'd', '*' '+' 'x' 's'} ;
% % colors = {'r',  [0 0.7 0], 'b', 'm', [0 0.7 0.7], [0.7 0 0.7], [0.5 0.5 0.7]};
% % time = (DT:DT:t)/day;
% optp = { [2 2 3 4 5 6 7], [1 2], 5, [1 2 3 4], [1.8] };
% Plot_val(preport.extras.DA,'savename','Eigenvalues','optp',optp,'x_lab',...
%     'Eigs','y_lab', 'log(eig)', 'revx', true,'o_legend', [{'kr_1', 'kr_2'}],...
%      'dir', [], 'figure', 10, 'titl', 'Relperm','o_ax', 2, 'W', W);
% 


linestyles = {'none', '-',  '--', '-.', ':' };
markerstyles = {'none', '.', 'o', 'd', '*' '+' 'x' 's'} ;
colors = {'r',  [0 0.7 0], 'b', 'm', [0 0.7 0.7], [0.7 0 0.7], [0.5 0.5 0.7]};
% optp is a vector containing the options to plot
% the elements are arranged as follows:
% 1. Marker, 2. Color, 3. Markersize, 4. Linestyle, 5. LineWidth
% Example
% optp = { [1 2 3 4 5 6 7], [1 2], 5, [1 2 3 4], [1.8] };
% x = 1:100;
% y = sin(2*x);
% figure
% plot(x,y,'Marker',markerstyles{optp{1}(2)},'MarkerEdgeColor', ...
%     colors{optp{2}(2)},'MarkerSize',optp{3} ,  'Linestyle', ...
%      linestyles{optp{4}(1)}, 'linewidth', optp{5})


    %%
opt = struct('optp',  {{ [1 2 3 4 5 6 7], [1 2], 5, [1 2 3 4], [1.8] }}, ...
    'figure',      0, ...
    'savename', [], ...
    'W',  [], ...
    'titl',   [], ...
    'x_lab',  [], ...
    'revx', false, ...
    'y_lab',  [], ...
    'o_ax',  2, ...
    'x', [], ...
    'o_legend','', ...
    'dir', []);

[opt] = merge_options(opt, varargin{:});
optp         = opt.optp;
nf           = opt.figure;
savename      = opt.savename;
W            = opt.W;
titl         = opt.titl;
x_lab        = opt.x_lab;
revx         = opt.revx;
y_lab        = opt.y_lab;
o_ax         = opt.o_ax;
x            = opt.x;
o_legend     = opt.o_legend;
dir          = opt.dir;
ly = size(y,2);
figure(nf);
clf
file = [savename];
if (~isempty(x))    
    for i = 1 : ly
f = plot(x, y(:,i),'Marker',markerstyles{optp{1}(i)},'MarkerEdgeColor', ...
    colors{optp{2}(i)},'MarkerSize',optp{3} ,  'Linestyle', ...
     linestyles{optp{4}(1)}, 'linewidth', optp{5});
 hold on
    end
else     
    for i = 1 : ly
f = plot(y(:,i),'Marker',markerstyles{optp{1}(i)},'MarkerEdgeColor', ...
    colors{optp{2}(i)},'MarkerSize',optp{3} ,  'Linestyle', ...
     linestyles{optp{4}(1)}, 'linewidth', optp{5});   
 hold on
    end
end
if (revx)
set (gca,'xdir','reverse')
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
    else if o_ax == 3
     axis  equal tight
        else
           axis off equal tight 
        end
    end
end

if (~isempty(dir))
        savefigures(gcf, file, dir)
end