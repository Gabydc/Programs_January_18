function subplotcbspei(nf,clim,ts,Np,G,nz,time,x)
figure(nf); hold on; axis off;

T = 1;
% xp = [0.05 0.25 0.45 0.65 0.85];
% xw = [0.02 0.02 0.02 0.02 0.1]; %yp
% yp = [0.15 0.15 0.15 0.15 0.05]; %xw
% yw = [0.95 0.95 0.95 0.95 0.75];  %yl

xp = [0.13 0.45 0.13 0.45 0.75];
xw = [0.6 0.6 0.1 0.1 0.1];      %yp
yp = [0.34 0.34 0.34 0.34 0.05]; %xw
yw = [0.34 0.34 0.34 0.34 0.75];  %yl

% set(h,'position',[0.85 0.25 0.05 0.5])

for i = 1 : 4 
 ax(i)=axes('position', [xp(i) xw(i) yp(i) yw(i)]);
%ax(i)=subplot(2,2,i);
%daspect([1 2 1])
axis equal tight off
plotCellData(G, x(:,T),'LineStyle','none');
%if nz > 1; view(-25,20); end
heading = [sprintf('T = %d days',time(T))];
title([ heading])
 T = ceil(ts/Np);
end
h=colorbar;
set(h, 'Position', [xp(5) xw(5) yp(5) yw(5)])



for i=1:4
      pos=get(ax(i), 'Position');
      set(ax(i), 'Position', [pos(1) 0.87*pos(2) 0.8*pos(3) pos(4)]);
end
set(gca,'LooseInset',get(gca,'TightInset'))

end

% ax(1)=subplot(4,1,1);
% 
% imagesc(magic(5),[1 64])
% 
% ax(2)=subplot(4,1,2);
% 
% imagesc(magic(6),[1 64])
% 
% ax(3)=subplot(4,1,3);
% 
% imagesc(magic(7),[1 64])
% 
% ax(4)=subplot(4,1,4);
% 
% imagesc(magic(8),[1 64])
% 
% h=colorbar;
% 
% set(h, 'Position', [.8314 .11 .0581 .8150])
% 
% for i=1:4
% 
%       pos=get(ax(i), 'Position');
% 
%       set(ax(i), 'Position', [pos(1) pos(2) 0.85*pos(3) pos(4)]);
% 
% end


% ax(2)=axes('position',[0.25 0.02 0.15 0.95]);
% T = ceil(ts/Np);
% axis equal tight off
% heading = [sprintf('T = %d days',ceil(time(T)))];
% plotCellData(G, x(:,T),'LineStyle','none');
% if nz > 1; view(-25,20); end
% title([ heading])
% ax(3)=axes('position',[0.45 0.02 0.15 0.95]);
% T = 2*ceil(ts/Np);
% axis equal tight off
% heading = [sprintf('T = %d days',ceil(time(T)))];
% plotCellData(G, x(:,T),'LineStyle','none');
% if nz > 1; view(-25,20); end
% title([ heading])
% ax(4)=axes('position',[0.65 0.02 0.15 0.95]);
% T = ts;
% axis equal tight off
% heading = [sprintf('T = %d days',ceil(time(T)))];
% plotCellData(G, x(:,T),'LineStyle','none');
% if nz > 1; view(-25,20); end
% title([ heading])
% caxis(clim);
% set(ax,'CLim',clim);
% h=colorbar;
% set(h,'position',[0.85 0.25 0.05 0.5])

%set(h,'position',[0.85 0.15 0.05 0.68])