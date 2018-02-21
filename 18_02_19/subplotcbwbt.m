function subplotcbwbt(nf,clim,wbt,Np,G,nz,time,x)
figure(nf); hold on; axis off;
N=size(x,2);
T = 1;
ax(1)=axes('position',[0.05 0.55 0.35 0.40]);
axis tight off
heading = [sprintf('T = %d PVI',time(T))];
plotCellData(G, x(:,1),'LineStyle','none');
caxis(clim);
if nz > 1; view(-25,20); end
title([ heading])
ax(3)=axes('position',[0.45 0.55 0.35 0.40]);
T = (wbt/Np);
axis tight off
heading = [sprintf('T = %.2f PVI',T)];
plotCellData(G, x(:,ceil(N/Np)),'LineStyle','none');
caxis(clim);
if nz > 1; view(-25,20); end
title([ heading])
ax(2)=axes('position',[0.05 0.05 0.35 0.40]);
T = 2*(wbt/Np);
axis tight off
heading = [sprintf('T = %.2f PVI',T)];
plotCellData(G, x(:,2*ceil(N/Np)),'LineStyle','none');
caxis(clim);
if nz > 1; view(-25,20); end
title([ heading])
ax(4)=axes('position',[0.45 0.05 0.35 0.40]);
T = 3*(wbt/Np);
axis tight off
heading = [sprintf('T = %.2f PVI',T)];
plotCellData(G, x(:,N),'LineStyle','none');
caxis(clim);
if nz > 1; view(-25,20); end
title([ heading])
caxis(clim);
set(ax,'CLim',clim);
h=colorbar;
set(h,'position',[0.85 0.05 0.08 0.90])