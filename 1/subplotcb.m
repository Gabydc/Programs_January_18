function subplotcb(nf,clim,ts,Np,G,nz,time,x)
figure(nf); hold on; axis off;

T = 1;
ax(1)=axes('position',[0.05 0.55 0.35 0.40]);
axis tight off
heading = [sprintf('T = %d days',time(T))];
plotCellData(G, x(:,T),'LineStyle','none');
if nz > 1; view(-25,20); end
title([ heading])
ax(2)=axes('position',[0.05 0.05 0.35 0.40]);
T = round(ts/Np);
axis tight off
heading = [sprintf('T = %d days',ceil(time(T)))];
plotCellData(G, x(:,T),'LineStyle','none');
if nz > 1; view(-25,20); end
title([ heading])
ax(3)=axes('position',[0.45 0.55 0.35 0.40]);
T = 2*round(ts/Np);
axis tight off
heading = [sprintf('T = %d days',ceil(time(T)))];
plotCellData(G, x(:,T),'LineStyle','none');
if nz > 1; view(-25,20); end
title([ heading])
ax(4)=axes('position',[0.45 0.05 0.35 0.40]);
T = 3*round(ts/Np)-1;
axis tight off
heading = [sprintf('T = %d days',ceil(time(T)))];
plotCellData(G, x(:,T),'LineStyle','none');
if nz > 1; view(-25,20); end
title([ heading])
caxis(clim);
set(ax,'CLim',clim);
h=colorbar;
set(h,'position',[0.85 0.05 0.08 0.90])