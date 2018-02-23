%...  The Matmol group (2016)
%...
    r = 0:0.01:3;
    y = min(2*r,2);
    area(r,y,'HandleVisibility','off')
    colormap spring
    ylabel('\phi(r)')
    xlabel('r')
    hold on
%...
    x = [0 1 3 3];
    y = [0 1 1 0];
    fill(x,y,'y','HandleVisibility','off');
    x2 = [0.5 1 2 1];
    y2 = [1   2 2 1];
    fill(x2,y2,'y','HandleVisibility','off')
%...
    for i=1:length(r),
        phi_van_leer(i)=(r(i)+abs(r(i)))/(1+abs(r(i)));
        phi_mc(i)=max(0,min([2*r(i) (1+r(i))/2 2]));
        phi_min_mod(i)=max(0,min([1 r(i)]));
        phi_superbee(i)=max([0 min(2*r(i),1) min(r(i),2)]);
        phi_smart(i)=max(0,min([4*r(i) (1+3*r(i))/4 2]));
        phi_koren(i)=max(0,min([2*r(i) (1+2*r(i))/3 2]));
    end
    plot(r(1:20:length(r)),phi_van_leer(1:20:length(r)),'--k')
    plot(r(1:20:length(r)),phi_mc(1:20:length(r)),'-.k')
    plot(r(1:20:length(r)),phi_min_mod(1:20:length(r)),'^k')
    plot(r(1:20:length(r)),phi_superbee(1:20:length(r)),'ok')
    plot(r(1:20:length(r)),phi_smart(1:20:length(r)),'*k')
    plot(r(1:20:length(r)),phi_koren(1:20:length(r)),'+k')
    legend('2nd-order TVD','TVD','','Van Leer','MC','Minmod','Superbee','Smart','Koren')
    xlabel('r')
    ylabel('\phi(r)')
    axis([0 3 0 2.5])

