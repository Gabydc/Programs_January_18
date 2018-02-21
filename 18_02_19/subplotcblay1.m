function subplotcblay1(nf,clim,ts,Np,G,nz,time,x)
figure(nf); hold on; axis off;
px = [0.05 0.25 0.45 0.65];
T = 1;
subplot(1,4,1,'Position',[px(1) 0.02 0.2 0.95]);
heading = [sprintf('T = %d days',time(T))]; 
title([ heading])
plotCellData(G, x(:,T),'LineStyle','none');
if nz > 1; view(3); daspect([ 6 3 0.5]); end
daspect([ 6 3 0.5]);
axis off 
for plotNo = 1 : Np-1
 T = (plotNo+1)*ceil(ts/(Np));
subplot(1,4,plotNo+1,'Position',[px(plotNo+1) 0.02 0.2 0.95]);
heading = [sprintf('T = %d days',time(T))]; 
title([ heading])
plotCellData(G, x(:,T),'LineStyle','none');
if nz > 1; view(3); daspect([ 6 3 0.5]); end
daspect([ 6 3 0.5]);
axis off 
 
end

h=colorbar;
set(h,'position',[0.85 0.13 0.05 0.69])

