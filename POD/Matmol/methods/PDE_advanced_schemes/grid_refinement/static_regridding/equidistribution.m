%...  The Matmol group (2016)
%...
    close all
    clear all
    
    global s
    z0 = -30.0;
    zL = 30.0;
    nz = 1001;
    dz = (zL-z0)/(nz-1);
    z  = [z0:dz:zL]';
    s  = 0.5;
    x  = kdv3_exact(z,0);
    
    plot(z,x)
    hold on
    nz2 = 21;
    dz2 = (zL-z0)/(nz2-1);
    z2  = [z0:dz2:zL]';
    x2  = kdv3_exact(z2,0);
    stem(z2,x2,'color','r')
    axis([z0 zL 0 0.3])
    xlabel('z')
    ylabel('x(z,t)')
    
    %...
    [arclen,seglen] = arclength(z2,100*x2);
    delta = arclen/nz2;
    mon(1)=0;
    for i = 1:length(seglen)
        mon(i+1) = mon(i)+seglen(i);
    end
    [z3,nz3] = invint(z2,mon,delta);
    x3 = kdv3_exact(z3,0);
    
    figure(2)
    plot(z,x)
    hold on
    stem(z3,x3,'color','r')
    axis([z0 zL 0 0.3])
    xlabel('z')
    ylabel('x(z,t)')
    
    [arclen,seglen] = arclength(z3,100*x3);
