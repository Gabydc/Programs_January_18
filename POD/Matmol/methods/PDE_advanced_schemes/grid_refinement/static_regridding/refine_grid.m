%... The Matmol group (2016)
    function status = refine_grid(t,x,flag,varargin)
%... this function refines the spatial grid after each time step,
%... interpolates the current values of the dependent variables,
%... and computes new differentiation matrices
%...
     global mu
     global z0 zL z nz D1 D2;
     global npdes nzmax alpha beta tolz bound imesh ilim
%...
%... refine the grid
%...
     [z_new,nz_new,ier,tolz_new]=agereg(z,x,npdes,nzmax,alpha,beta,tolz,bound,imesh,ilim);
%...
     figure(1)
     hold on
     plot(t,nz_new)
%...
%... interpolate the dependent variables
%...
     z=z_new';
     nz=nz_new;
     tolz=tolz_new;
%...
     x_new=spline(z,x,z_new);
     x=x_new';
%...
     figure(2)
     subplot(2,1,1)
     plot(z,x);
     subplot(2,1,2)
     plot(z,t*ones(nz,1),'.')
%...
%... compute new differentiation matrices
%...
%     D1=budm1d1nu(z,v);
%     D1=budm1d3nu(z,v);
     D1=budm1d4nu(z,v);
%...
     D2=cdm2d2nu(z);
%     D2=cdm2d4nu(z);
%...
