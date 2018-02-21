for i = 1 : 10
    cnn(i) = preport(i).cond;
end
for i = 11 : ts
    cnn(i) = preport(1,i).cond(1);
end


    f(5) = figure(5);
    if def == 0
     
            hn=plot((dT:dT:T)/day,cnn(1:ts),'r*');
               axis square
    ylabel(' \kappa_2(M^{-1}A)','FontSize',16)
    xlabel('Time (days) ','FontSize',16)
        axis square
    file{5} = 'condest';
   B=[dir2   file{5}   '.jpg'];
                 saveas(f(5),B) 
    else
        
    subplot(1,2,1)
            hn=plot((dT:dT:10*dT)/day,cnn(1:10),'r*');
               axis square
    ylabel(' \kappa_2(M^{-1}A)','FontSize',16)
    xlabel('Time (days) ','FontSize',16)
            subplot(1,2,2)
            hn=plot((11*dT:dT:T)/day,cnn(11:ts),'b*');
       % legend('ICCG');
        axis square
    ylabel(' \kappa_2(PM^{-1}A)','FontSize',16)
    xlabel('Time (days) ','FontSize',16)
    file{5} = 'condest';
   B=[dir2   file{5}   '.jpg'];
                 saveas(f(5),B)
          f(6) = figure(6);
file{6} = 'condestdef';
            hn=plot((dT:dT:10*dT)/day,cnn(1:10),'r*');
            hold on    
            hn=plot((11*dT:dT:T)/day,cnn(11:ts),'b*');
       % legend('ICCG');
        axis square
    ylabel(' \kappa_2','FontSize',16)
    xlabel('Time (days) ','FontSize',16)
    file{6} = 'condest_1';
   B=[dir2   file{6}   '.jpg'];
                 saveas(f(6),B)
          
          
    end