%... The MatMol Group (2016)
%... Evolution of Ak in the matrix method. Von Neumann stability
%... analysis 

    close all
    clear all
    
%... Values of sigma    
    sigma=[.8 1.2 1.9 2.1];
    n=20;
    puis = [1:1:99 100:10:990 1000:100:5900];
    
%... Compute Ak and plot the results    
    figure
    hold
    for i=1:length(sigma)
        a=diag((1-sigma(i))*ones(n,1),0)+diag(sigma(i)*ones(n-1,1),-1);
        for k=1:length(puis)
            ak(k,i)= norm(a^puis(k));
        end
    end
    plot(log10(puis),log10(ak(:,1)),'color',[1 0 0],'linewidth',2)
    plot(log10(puis),log10(ak(:,2)),'color',[0 1 0],'linewidth',2)
    plot(log10(puis),log10(ak(:,3)),'color',[0 0 1],'linewidth',2)
    plot(log10(puis),log10(ak(:,4)),'color',[.5 0 0.5],'linewidth',2)

    legend('\sigma = 0.8','\sigma = 1.2','\sigma = 1.9','\sigma = 2.1','Location','Northwest')
    title('||A^k|| = f(k,\sigma)')

    set(gca,'Ytick',-40:10:40)
    set(gca,'YTickLabel',{'1e-40','1e-30','1e-20','1e-10','1e+0','1e+10','1e+20','1e+30','1e+40'})
    set(gca,'XTick',0:1:2)
    set(gca,'XTickLabel',{'1e+0','1e+1','1e+2'})
    axis([-.1 2.1 -42 42])
    grid on
    