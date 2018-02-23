%...  The MatMol Group (2016)
    function plot_solar(y)

% plot solar system

    figure('units','normalized','outerposition',[0 0 1 1])
    hold on
    grid on
    Ycorr = [y(1,:); y(2,:); y(3,:); y(1,:); y(2,:); y(3,:);...
        y(1,:); y(2,:); y(3,:); y(1,:); y(2,:); y(3,:);...
        y(1,:); y(2,:); y(3,:); y(1,:); y(2,:); y(3,:)];
    y(1:18,:) = y(1:18,:) - Ycorr;
    plot3(y(1,:),y(2,:),y(3,:),'g')     % sun
    plot3(y(4,:),y(5,:),y(6,:),'m')     % Jupiter
    plot3(y(7,:),y(8,:),y(9,:),'r')     % Saturn
    plot3(y(10,:),y(11,:),y(12,:),'c')  % Uranus
    plot3(y(13,:),y(14,:),y(15,:),'b')  % Neptune
    plot3(y(16,:),y(17,:),y(18,:),'k')  % Pluto
    
    plot3(y(1,1),y(2,1),y(3,1),'g.','MarkerSize',20);
    plot3(y(4,1),y(5,1),y(6,1),'m.','MarkerSize',20);
    plot3(y(7,1),y(8,1),y(9,1),'r.','MarkerSize',20);
    plot3(y(10,1),y(11,1),y(12,1),'c.','MarkerSize',20);
    plot3(y(13,1),y(14,1),y(15,1),'b.','MarkerSize',20);
    plot3(y(16,1),y(17,1),y(18,1),'k.','MarkerSize',20);
    
    text(y(1,1),y(1,1)+1,y(3,1),'S','Fontsize',18);
    text(y(4,1),y(5,1)+1,y(6,1)+0.5,'J','Fontsize',18);
    text(y(7,1)+0.5,y(8,1)+1,y(9,1),'S','Fontsize',18);
    text(y(10,1),y(11,1)+1.2,y(12,1),'U','Fontsize',18);
    text(y(13,1),y(14,1)+1.2,y(15,1),'N','Fontsize',18);
    text(y(16,1)+0.3,y(17,1)+1,y(18,1),'P','Fontsize',18);
    
    set(gca,'CameraPosition',[-214.0535 -786.1018  -81.0481])
    