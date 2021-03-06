nf = 0;

pmark =['o' '+' '*' 'x' 'd' 's' ];
pcol ={'r' [0 0.7 0] 'b'  [0 0.7 0.7] [0.7 0 0.7] [0.5 0.5 0.7] };
time = (DT:DT:t)/day;

% Plot permeability
nf = nf + 1;
figure(nf);
file{nf} = ['Permeability'];
perm = rock.perm(:,1);
permD = log10(perm/darcy);
plotCellData(G, permD,'LineStyle','none');
if nz>1; view(-40,20);  daspect([ 3 3 3]);  end
axis off
daspect([ 3 3 3]);
if(use_wells)
    plotWell(G, W(1:4),'color', 'r')
    plotWell(G, W(5), 'color', 'b')
end
caxis([min(permD) max(permD)])
c = colorbar('southoutside');
c.Label.String = '\fontsize{12} log(Permeability) [D]';
%%
if(use_wells)
    plotWell(G, W(1:4),'color', 'r')
    plotWell(G, W(5), 'color', 'b')
    
    % Plot Pressures
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
    
    nf = nf + 1;
    figure(nf);
    file{nf} = ['Pressure_wp'];
    for i=1:numel(W)-1
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

%% Pressure
pmin=0.9*min(Pressure(:,2)/barsa);
pmax=1.1*max(Pressure(:,k)/barsa);

nf = nf + 1;
figure(nf);
file{nf} = ['Pressure'];
clim = [pmin pmax];
plotCellData(G, Pressure(:,k)/barsa,'LineStyle','none');
if(use_wells)
    plotWell(G, W(1:4),'color', 'r');
    plotWell(G, W(5), 'color', 'b');
end
if nz>1; view(3);  daspect([ 6 3 0.5]); end;
daspect([ 6 3 0.5]);
axis off
caxis(clim)
c = colorbar;
%c = colorbar('southoutside');
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
if nz>1; view(3);  daspect([ 6 3 0.5]);  end
axis off
daspect([ 6 3 0.5]);
caxis([0 1])
%c = colorbar('southoutside');
c = colorbar;
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
if nz>1; view(3);  daspect([ 6 3 0.5]);  end
axis off
daspect([ 6 3 0.5]);
caxis([0.5 1])
%c = colorbar('southoutside');
c = colorbar;
c.Label.String = '\fontsize{12} Saturation';






%%
if use_DICCG
    nf = nf+1;
    file{nf} = ['eig_pod'];
    f(nf) = figure(nf);
    plot(log((diag(S))),'*r');
    set(gca, 'XDir','reverse')
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
% Plot the number of iterations
for j = 1 : k
    its(j)=preport(j).iter;
end
nf = nf + 1;
f(nf) = figure(nf);
file{nf} = ['Iterations'];
if(use_ICCG)
    plot(time,its(1:k),'r*');
    legend('ICCG');
else
    if(window)
        plot(time(1:length(dpod)),its(1:length(dpod)),'r*');
        hold on
        plot(time(length(dpod)+1:k),its(length(dpod)+1:k),'b*');
        legend('ICCG','DICCG');
    else
        if(training)
            plot(time,its(1:k),'r*');
            legend('ICCG');
        else
            plot(time,its(1:k),'b*');
            legend('DICCG');
        end
        
    end
end


axis([ DT/day t/day 0 max(its)]);
ylabel('Number of iterations','FontSize',16)
xlabel('Time (days)','FontSize',16)
axis square

%%

Np = 4;
nf = nf + 1;
figure(nf);
file{nf} = ['Water_saturation'];
clim = [0 1];
if(model_SPE)
    subplotcbspe(nf,clim,k,Np,G,nz,time,Sat)
else
    subplotcbspe(nf,clim,k,Np,G,nz,time,Sat)
end
%%

nf = nf + 1;
figure(nf);
file{nf} = ['4_Pressure'];
clim = [pmin pmax];
subplotcbspe(nf,clim,k,Np,G,nz,time,Pressure/barsa)




for i = 1 : nf
    f(i) = figure(i);
    %set(f(i), 'Visible', 'off')
    savefigures(f(i), file{i}, dir2)
end
clear figure
clear f




