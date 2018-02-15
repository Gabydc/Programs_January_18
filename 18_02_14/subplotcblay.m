function subplotcblay(nf,clim,ts,Np,G,nz,time,x)
figure(nf); hold on; axis off;

T = 1;
for
h=subplot(1,4,plotNo+1,'Position',[px(plotNo+1),0.03,0.17,0.8]);
ax(1)=axes('position',[0.05 0.02 0.15 0.95]);
axis equal tight off
heading = [sprintf('T = %d days',time(T))];
plotCellData(G, x(:,T),'LineStyle','none');
hb = colorbar('Position',[px(plotNo+1)+0.18,0.1,0.02,0.5]);



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
set(h,'position',[0.85 0.40 0.05 0.20])

%set(h,'position',[0.85 0.15 0.05 0.68])