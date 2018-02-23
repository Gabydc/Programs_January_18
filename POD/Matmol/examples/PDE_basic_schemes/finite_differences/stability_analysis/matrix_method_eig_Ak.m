%... The MatMol Group (2016)
%... Eigenvalues loci for different values of sigma in the Von
%... Neumann stability analysis

    close all
    clear all

%... Set sigma and the ange
    sigma = 0.2:0.2:1.4;
    th    = [0:0.001:2*pi 2*pi]';

%... Compute data
    for ii = 1 : length(sigma)
        vp = 1 - sigma(ii) + sigma(ii)*exp(j*th);
        re_vp(:,ii) = real(vp);
        im_vp(:,ii) = imag(vp);
    end

%... Create the figure
    figure
    hold
    plot(cos(th),sin(th),'linewidth',2)
    for i=1:length(sigma)
        vp=1-sigma(i)+sigma(i)*exp(j*th);
        plot(real(vp),imag(vp))
        text(.8-sigma(i)*2,sigma(1)/2,sprintf('%...3.1f',sigma(i)))
    end
    text(-2.3,sigma(1)/2,'\sigma = ')
    plot([-2.4 1.1],[0 0])
    plot([0 0],[-1.75 1.75])

    axis([-2.4 1.1 -1.75 1.75])
    axis square

