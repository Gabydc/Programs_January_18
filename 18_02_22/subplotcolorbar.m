function subplotcolorbar(clim,ts,Np,G,nz,time,x)
%clim = [pmin pmax];

figure; hold on; axis off;
T = 1;
heading = [sprintf('T=%.2d days',time(T))];
ax(1)=axes('position',[0.05 0.55 0.35 0.40]);
T = ceil(ts/Np);
axis tight off
heading = [sprintf('T=%.2d days',time(T))];
plotCellData(G, x(:,T),'LineStyle','none');
title([ heading])
ax(2)=axes('position',[0.05 0.05 0.35 0.40]);
T = 2*ceil(ts/Np);
axis tight off
heading = [sprintf('T=%.2d days',time(T))];
plotCellData(G, x(:,T),'LineStyle','none');
title([ heading])
ax(3)=axes('position',[0.45 0.55 0.35 0.40]);
T = 3*ceil(ts/Np);
axis tight off
heading = [sprintf('T=%.2d days',time(T))];
plotCellData(G, x(:,T),'LineStyle','none');
title([ heading])
ax(4)=axes('position',[0.45 0.05 0.35 0.40]);
T = ts;
axis tight off
heading = [sprintf('T=%.2d days',time(T))];
plotCellData(G, x(:,T),'LineStyle','none');
title([ heading])
caxis(clim);
set(ax,'CLim',clim);
h=colorbar;
set(h,'position',[0.85 0.05 0.10 0.90])