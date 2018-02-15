function subplotcbspe(nf,clim,ts,Np,G,nz,time,x)
figure(nf); hold on; axis off;

T = 1;
xp = [0.05 0.25 0.45 0.65 0.95];
xw = [0.02 0.02 0.02 0.02 0.25];
yp = [0.15 0.15 0.15 0.15 0.05];
yw = [0.95 0.95 0.95 0.95 0.5];
for i = 1 : 4
ax(i)=axes('position',[0.05 0.02 0.15 0.95]);
axis equal tight off
heading = [sprintf('T = %d days',time(T))];
plotCellData(G, x(:,T),'LineStyle','none');
if nz > 1; view(-25,20); end
title([ heading])
ax(2)=axes('position',[0.25 0.02 0.15 0.95]);
T = ceil(ts/Np);
axis equal tight off
heading = [sprintf('T = %d days',ceil(time(T)))];
plotCellData(G, x(:,T),'LineStyle','none');
if nz > 1; view(-25,20); end
title([ heading])
ax(3)=axes('position',[0.45 0.02 0.15 0.95]);
T = 2*ceil(ts/Np);
axis equal tight off
heading = [sprintf('T = %d days',ceil(time(T)))];
plotCellData(G, x(:,T),'LineStyle','none');
if nz > 1; view(-25,20); end
title([ heading])
ax(4)=axes('position',[0.65 0.02 0.15 0.95]);
T = ts;
axis equal tight off
heading = [sprintf('T = %d days',ceil(time(T)))];
plotCellData(G, x(:,T),'LineStyle','none');
if nz > 1; view(-25,20); end
title([ heading])
caxis(clim);
set(ax,'CLim',clim);
h=colorbar;
set(h,'position',[0.85 0.25 0.05 0.5])

%set(h,'position',[0.85 0.15 0.05 0.68])