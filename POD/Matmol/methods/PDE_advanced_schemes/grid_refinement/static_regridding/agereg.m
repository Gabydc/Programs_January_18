%... The Matmol group (2016)
    function [z_new,nz_new,ier,tolz_new]=agereg(z,x,npdes,nzmax,alpha,beta,tolz,bound,imesh,ilim)
%... agereg computes a new grid, which is subequidistributing a monitor function 
%...                             and which is locally bounded
%...
%... This implementation is based on AGE (a previous code by P. Saucez, W.E. Schiesser 
%... and A. Vande Wouwer - Journal of Computational Physics, 1996) and an algorithm 
%... by Kautsky and Nichols for equidistribution with constraints (). AGEREG builds upon
%... NEWMESH, a code developed by Gerd Steinebach (Ph.D. Thesis, 1993)
%... Modifications by P. Saucez, W.E. Schiesser and A. Vande Wouwer (2001)  
%...
%... input parameters:
%...
%...     z(nz)          : current location of nodes
%...     x(npdes*nz)    : current dependent variable matrix
%...     npdes          : number of PDEs
%...     nzmax          : maximum number of nodes
%...     alpha          : parameter of the monitor function 
%...                     (limit the maximum grid interval, i.e., 
%...                      hmax**2 < tolz**2/alpha)
%...     beta           : parameter of the monitor function
%...                      (limit the excessive clustering of nodes) 
%...     tolz           : reference quantity for the equidistribution 
%...                      of the monitor function
%...     bound >=1      : parameter defining a locally bounded grid
%...     imesh      =0  : monitorfunction f = dqrt( alpha + sum (x_z)^2 )
%...                =1  : monitorfunction f = dsqrt( alpha + max(dabs(x_zz)))
%...     ilim       =0  : no limit on the second spatial derivatives of x
%...                =1  : limit the value of x_zz to beta 
%...
%... output parameters:
%...
%...     nz_new         : new number of nodes
%...     z_new(nz_new)  : new location of nodes
%...     ier      = 0   : no error
%...              = 1   : warning, tolz is too small
%...     tolz_new       : adjusted value of tolz
%...
%***************************************************************************
%...
%...  initialisation
      nz=length(x)/npdes;
      tolz_new=tolz;
      ier=0;
%...
%...  Definition of the monitor function
%...
%...  compute in xz and xzz the value of the first and second 
%...  spatial derivatives of each dependent variable using cubic 
%...  splines
      for i=1:npdes;      
          xval(1:nz)=x(i:npdes:npdes*nz-(npdes-i));
          s=ncs(z',xval);
          xz(i,1:nz)=s(1:nz,3)';
          xzz(i,1:nz)=2*s(1:nz,2)';             
      end;          
%...
%...	depending on imesh, the monitor function is in the form
%...	sqrt(alpha+sum(x_z^2)) or sqrt(alpha+max(abs(x_zz)))
%...
%...  if imesh = 1, determine max(abs(x_zz))
	   if imesh == 1
          if npdes == 1
          xzzm=abs(xzz);
          else
          xzzm=max(abs(xzz));
          end
%...  
%...      if ilim = 1, limit the value of max(abs(x_zz)) 
	       if ilim == 1
                 index = (xzzm > beta*ones(1,nz));
                 xzzm(index) = beta;
          end
      end     
%...
%...  assemble the monitor function
%...
      if imesh == 0
          mon=sum(xz.^2,1);
      elseif imesh == 1
	      mon=xzzm;
	  end
      fmon=sqrt(alpha+mon);
%...
%...  evaluation of the integral of the monitor function using
%...  the trapezoidal rule
%...
      [fmoni]=trapez(z,fmon);
%...
%...  determination of delta (which is the quantity with respect 
%...  to which the padded function is to be equidistributed) taking
%...  into account the maximum number of nodes
%...       
%...  first, consider the integral of the monitor function and
%...  determine the number of nodes required for an equidistribution 
%...  with respect to tolz_new=tolz
%...
      nz_new=fmoni(nz)/tolz_new;
%...
         if nz_new > (nzmax-1)
%...
%...     tolz is too small and should be increased accordingly
%...     in turn, nzmax nodes are used
%...
            tolz_new=fmoni(nz)/(nzmax-1);
	        nz_new=nzmax;
            ier=1
%...       
%...  otherwise round the number of nodes
%...
         else
            nz_new=ceil(nz_new);
         end
%...
%...  there must be at least 10 nodes
%...
      nz_new=max(nz_new,10);
%...
%...  estimation of delta and lambda
%...
      delta=fmoni(nz)/(nz_new-1);
      lambda=log(bound)/delta;
%...
%...  padd the monitor function 
%...
      [fmon]=padding(z,fmon,lambda);
%...
%...  compute the integral of the padding function using 
%...  the trapezoidal rule
%...
      [fmoni]=trapez(z,fmon);
%...
%...  second, consider the integral of the padding function and
%...  determine the number of nodes required for an equidistribution 
%...  with respect to tolz_new
%...
      nz_new=fmoni(nz)/tolz_new;
%...
          if nz_new > (nzmax-1)
%...
%...      tolz_new is too small and should be increased accordingly
%...      in turn, nzmax grid points are used
%...
             tolz_new=fmoni(nz)/(nzmax-1);
	         nz_new=nzmax;
             ier=1
%...       
%...  otherwise round the number of nodes
%...
          else
             nz_new=ceil(nz_new);
          end
%...
%...  there must be at least 10 nodes
%...
      nz_new=max(nz_new,10);
%...
%...  second estimation of delta
%...
      delta=fmoni(nz)/(nz_new-1);
%...
%...  check the resulting value of bound
%...
      bound1=exp(delta*lambda)
%...
%...
%...  determination of the new mesh by inverse linear interpolation
%...
      [z_new,nz_new]=invint(z,fmoni,delta);

