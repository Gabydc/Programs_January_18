nf = 0;

pmark =['o' '+' '*' 'x' 'd' 's' ];
pcol ={'r' [0 0.7 0] 'b'  [0 0.7 0.7] [0.7 0 0.7] [0.5 0.5 0.7] };


% Plot permeability
nf = nf + 1;
figure(nf);
file{nf} = ['Permeability'];  
perm = rock.perm(:,1);
plotCellData(G, log(perm),'LineStyle','none');
 view(3)
 axis equal off
 daspect([ 3 3 1])
%c = colorbar();
%title ('\fontsize{12} Permeability field [mD]','southoutside');
caxis([min(log(perm(:))) max(log(perm(:)))])
if(use_wells)
 plotWell(G, W(1:4),'color', 'r')
 plotWell(G, W(5), 'color', 'b')

%% Plot Pressures
% Wells pressures
nf = nf + 1;
figure(nf);
file{nf} = ['Pressure_w'];   
for i=1:5
    for j = 1 : nstep
        if training
    wi(j,i)=W1{j}(i).val/barsa;
        else
            wi(j,i)=W0{j}(i).val/barsa;
        end
    end
    plot(wi(:,i),'color', pcol{i})
    hold on
end
xlabel('Time [d]'), ylabel('Pressure [bars]')
legend({ W(1:end).name})
%%
nf = nf + 1;
figure(nf);
file{nf} = ['Pressure_wp'];
for i=1:5
    plot(pw1(i,:)/barsa,'color', pcol{i},'Marker',pmark(i))
    hold on
end
xlabel('Time [d]'), ylabel('Pressure [bars]')
legend({ W(1:end-1).name })

nf = nf + 1;
figure(nf);
file{nf} = ['Tot_prod_rate'];
plot(convertTo(Prod.t, day), convertTo(Prod.vpt(:,1:end-1), meter^3/day))
legend({ W(1:end-1).name }, 'Location', 'Best')
xlabel('Time [d]'), ylabel('Total Production Rate [m^3/d]')



nf = nf + 1;
figure(nf);
file{nf} = ['Oil_prod_rate'];
plot(convertTo(Prod.t, day), convertTo(Prod.opr(:,1:end-1), meter^3/day))
legend({ W(1:end-1).name }, 'Location', 'Best')
xlabel('Time [d]'), ylabel('Oil Production Rate [m^3/d]')

nf = nf + 1;
figure(nf);
file{nf} = ['Water_prod_rate'];
plot(convertTo(Prod.t, day), convertTo(Prod.wpr(:,1:end-1), meter^3/day))
legend({ W(1:end-1).name }, 'Location', 'Best')
xlabel('Time [d]'), ylabel('Water Production Rate [m^3/d]')

nf = nf + 1;
figure(nf);
file{nf} = ['Well_wat_cut'];
plot(convertTo(Prod.t,day), Prod.wc(:,1:end-1))
legend({ W(1:end-1).name }, 'Location', 'Best')
xlabel('Time [d]'), ylabel('Well Water Cut')
end
% Solution
pmin=0.9*min(Pressure(:,2)/barsa);
pmax=1.1*max(Pressure(:,k)/barsa);

nf = nf + 1;
figure(nf);
%title('Pressure [bars]');
file{nf} = ['Pressure'];
clim = [pmin pmax];
plotCellData(G, Pressure(:,k)/barsa,'LineStyle','none');
if(use_wells)
plotWell(G, W(1:4),'color', 'r');
plotWell(G, W(5), 'color', 'b');
end
view(3)
axis equal off
daspect([ 3 3 1])
caxis(clim)
c = colorbar('southoutside');
c.Label.String = '\fontsize{12} Pressure [bars]';

%% Saturation



nf = nf + 1;
figure(nf);
file{nf} = ['Sat'];
plotCellData(G, x.s(:,1), 'EdgeColor', 'k', ...
    'EdgeAlpha', 0.050, 'FaceAlpha', 0.375,'LineStyle','none')
if(use_wells)
plotWell(G, W(1:4),'color', 'r');
plotWell(G, W(5), 'color', 'b');
end
view(3)
axis equal off
daspect([ 3 3 1])
caxis([0 1])
c = colorbar('southoutside');
c.Label.String = '\fontsize{12} Saturation';



nf = nf + 1;
figure(nf);
file{nf} = ['Sat_1'];
plotCellData(G, x.s(:,1), find(x.s > 0.5), 'EdgeColor', 'k', ...
    'EdgeAlpha', 0.050, 'FaceAlpha', 0.375,'LineStyle','none')
if(use_wells)
plotWell(G, W(1:4),'color', 'r');
plotWell(G, W(5), 'color', 'b');
end
view(3)
axis equal off
daspect([ 3 3 1])
caxis([0.5 1])
c = colorbar('southoutside');
c.Label.String = '\fontsize{12} Saturation';





%%

% Press = x.pressure;
% Sat_w = x.s;
% 
% 
% plotNo = 0;
% Np = 3;
% time = (0:DT:DT*k)/day;
% px = [0.01 0.26 0.51 0.76];
% 
% nf = nf + 1;
% figure(nf);
% %title('Water Saturation');
% file{nf} = ['Water_saturation'];
% clim = [0 0.75];
% nz = numel(layers);
% subplot(2,2,1,'Position',[0.05 .5 0.4 0.4])
% plotCellData(G, Sat_w,'LineStyle','none', 'EdgeColor', 'k', ...
%     'EdgeAlpha', 0.050, 'FaceAlpha', 0.375);
% % for i = 1 : 4
% %     plotWell(G, W(i),'color', 'r');
% % end
% % plotWell(G, W(5), 'color', 'b');
% view(3)
% axis equal off
% % colorbar('southoutside') 
% % caxis(clim)

%%
if use_DICCG
    nf = nf+1;
    file{nf} = ['eig_pod'];
    f(nf) = figure(nf);
    plot(log((diag(S))),'*r');
    ylabel('log(Value) ','FontSize',16)
    xlabel('Eigenvalue','FontSize',16)
    axis('tight');
    nf = nf + 1;
    
    file{nf} = ['eig_vect' ];
    f(nf) = figure(nf);
    plot(Z);
    ylabel('vectors ','FontSize',16)
    xlabel('eigenvectors','FontSize',16)
    axis('tight')
    end


%%



% subplotcbspe(nf,clim,k,Np,G,nz,time,Sat_w)
% subplotcbwbt(nf,clim,wbt,Np,G,nz,time,Sat_w)

% subplotcbspe(nf1,clim,k,Np,G,nz,time,Press)

for i = 1 : nf
    f(i) = figure(i);
    savefigures(f(i), file{i}, dir2)
end
clear figure
clear f




