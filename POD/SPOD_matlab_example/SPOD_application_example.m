%% SPOD application example

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Written by: Moritz Sieber
% Version: 03/2017
% E-Mail: moritz.sieber@fd.tu-berlin.de
%
% Please cite our article if you use this method in your own work:
%
% Sieber, M., Paschereit, C.O. and Oberleithner, K. (2016) ‘Spectral proper 
% orthogonal decomposition’, Journal of Fluid Mechanics, 792, pp. 798–828. 
% doi: 10.1017/jfm.2016.103.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear 

% load data
load('PIV_jext_example.mat')

% subtract mean
Um = mean(U,3);
Vm = mean(V,3);
U1 = bsxfun(@minus,U,Um);
V1 = bsxfun(@minus,V,Vm);

% plot mean magnitude
figure
contourf(X,Y,sqrt(Um.^2+Vm.^2),'edgecolor','none')
hold on
h = streamslice(double(X),double(Y),double(Um),double(Vm));
set(h,'color','k')

xlabel('x')
ylabel('y')
title('mean velocity field')
legend('velocity magnitude','streamlines')

%% calculate SPOD
Nfilt = 25;    % filter size
Npod = [];     % number of modes in output
Wxyz = [];     % spatial weighting
boundary = 1;  % periodic BC
corr_type = 1; % snapshot correlation
[Ui,Vi,a,lambda,norma] = SPOD(U1,V1,Nfilt,Npod,Wxyz,boundary,corr_type);

% find linked modes
[~,idx] = min(abs(cumsum(lambda)/sum(lambda)-0.95));
mode = SPODpost(a(:,1:idx));

% plot SPOD spectrum
h = figure;
pos = get(h,'position');
pos(1) = pos(1)+100; pos(2) = pos(2)-150;
set(h,'position',pos);

scale = 50;
f = [mode.f]; K = [mode.K]*100; c = [mode.c];
scatter(f,K,c*scale+0.01,c,'filled')
set(gca,'xscale','lin','yscale','log')
xlabel('frequency [1/sample]')
ylabel('energy [% of TKE]')
box on
for i=1:4
    text(f(i),K(i),num2str(i),...
        'VerticalAlignment','bottom')
end
title('SPOD spectrum')

%% plot SPOD modes
h = figure;
pos = get(h,'position');
pos(1) = pos(1)+200; pos(2) = pos(2)-300;
set(h,'position',pos);

for i=1:4
    subplot(4,3,(i-1)*3+1)
    idx = i;
    contourf(X,Y,Vi(:,:,mode(idx).ind(1)),10,'edgecolor','none')
    text(1,13,sprintf('\\Phi_{v,%i,1}',i))
    
    subplot(4,3,(i-1)*3+2)
    contourf(X,Y,Vi(:,:,mode(idx).ind(2)),10,'edgecolor','none')
    text(1,13,sprintf('\\Phi_{v,%i,2}',i))
    if i==1
        title('SPOD mode shapes and coefficient spectra')
    end
    
    subplot(4,3,(i-1)*3+3)
    psd = abs(fft(mode(idx).a)).^2;
    nf = round(length(psd)/2)+1;
    f = linspace(0,0.5,nf);
    semilogy(f(2:nf),psd(2:nf))
    axis([0 0.5 1e-15 1e2])
    text(0.42,0.5,sprintf('a_{%i}',i))
end
