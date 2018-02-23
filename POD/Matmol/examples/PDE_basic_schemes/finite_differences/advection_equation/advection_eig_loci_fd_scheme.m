%... The MatMol Group (2016)
%... Eigenvalue loci for the advection equation with different
%... finite difference schemes
    close all
    clear all

%... Angle
    th = [0:.001:2*pi 2*pi];

%... 3-point centered
    vp_center    = -j*100*sin(th);
    re_vp_center = real(vp_center);
    im_vp_center = imag(vp_center);

%... 2-point downwind
    vp_down    = 100*(1-cos(th)-j*sin(th));
    re_vp_down = real(vp_down);
    im_vp_down = imag(vp_down);

%... 2-point upwind
    vp_up    = 100*(cos(th)-1-j*sin(th));
    re_vp_up = real(vp_up);
    im_vp_up = imag(vp_up);

%... Plot the results 
%... 3-point centered
    figure
    hold
    plot(re_vp_center,im_vp_center,'linewidth',2)
    text(-200,200,'\lambda_k spectrum')
    text(-200,180,'3-point centered')
    plot([-220 220],[0 0])
    plot([0 0],[-220 220])
    axis([-220 220 -220 220])
    axis square

%... 2-point downwind    
    figure
    hold
    vp=100*(1-cos(th)-j*sin(th));
    plot(re_vp_down,im_vp_down,'linewidth',2)
    text(-200,200,'\lambda_k spectrum')
    text(-200,180,'2-point downwind')
    plot([-220 220],[0 0])
    plot([0 0],[-220 220])
    axis([-220 220 -220 220])
    axis square

%... 2-point upwind    
    figure
    hold
    vp=100*(cos(th)-1-j*sin(th));
    plot(re_vp_up,im_vp_up,'linewidth',2)
    text(-200,200,'\lambda_k spectrum')
    text(-200,180,'2-point upwind')
    plot([-220 220],[0 0])
    plot([0 0],[-220 220])
    axis([-220 220 -220 220])
    axis square
