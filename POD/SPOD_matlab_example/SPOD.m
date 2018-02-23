function varargout = SPOD(varargin)
% [Ui,Vi,Wi,...,a,lambda,mode_norm] = ...
%                       SPOD(U,V,W,...,Nfilt,Npod,Wxyz,boundary,corr_type)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [Ui,Vi,Wi,...,a,lambda,mode_norm] = ...
%                       SPOD(U,V,W,...,Nfilt,Npod,Wxyz,boundary,corr_type)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function calculates the Spectral Proper Orthogonal Decomposition.
% The SPOD provides basis vectors and Fourier coefficients of the input
% data.
% 
% The input data can be recomposed from the output data (U = Ui*a'):
% U(k,i) = Ui(k,1:Npod) * a(i,1:Npod)'
% If Npod < Nsnap the expression is only approximately equal
%
% Rquired function inputs: 
%
% U,V,W,... :   Arbitrary number of velocity components. Each array may 
%               have any number of dimensions but all must have the same
%               number of dimensions. The last dimension must span the
%               temporal direction of the data e.g. U(:,:,i) is the i-th
%               snapshot of a 2D data set.
%
% Optional function inputs:
% Use empty arguments to apply default values, e.g. SPOD(U,Nfilt,[],[],1)
%
% Nfilt:        The length of the SPOD filter. Scalar integer greater than  
%               or equal to 0. Default is 0 (classical POD). 
% Npod:         Number of modes to return. Scalar integer greater than 0.
%               Default is number of snapshots. 
% Wxyz:         Spatial weighting for the inner product. Double array with
%               the same dimensions as the spatial dimensions of U,V,W... .
%               The weighted inner product reads <u,v>_w = u'*diag(w)*v.
%               Default is uniform spatial weighting e.g.
%               Wxyz = ones(size(U(:,:,1)))
% boundary:     Handling of missing data at start and end of the time
%               series. Scalar with value:
%               0 -> pad with zeros, 
%               1 -> use periodic boundary conditions. 
%               Default is 1 (periodic).
% corr_type:    Type of correlation used to calculate the SPOD. Scalar with 
%               value:
%               0 -> spatial correlation (use when Ngrid*Nfilt < Nsnap), 
%               1 -> temporal correlation (use when Ngrid*Nfilt > Nsnap).
%               Default is 1(snapshot POD).
%               
% Function outputs:
%
% Ui,Vi,Wi,... :POD modes (spatial modes). Number of output modes depends
%               on the number input velocity components. The dimension of
%               the output data is the same as the input data, whereas the
%               last dimension contains individual modes instead of
%               snapshots.
% a:            POD coefficients (temporal modes). 2D array where the first
%               dimension represents temporal snapshots and the second the
%               single modes, e.g. a(:,1) is the time-series of the first
%               mode coefficient.
% lambda:       Energy of individual modes. 1D array of length Npod. 
%               diag(a'*a)/Nsnap = lambda.
% mode_norm:    Norm of the single spatial modes. 1D array of length Npod. 
%               For classical POD (Nfilt=0) the norm is 1 for all modes. 
%               sqrt(<Ui,Ui>_w + <Vi,Vi>_w + <Wi,Wi>_w + ...) = mode_norm. 
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
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
 
 
% parse input arguments %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[Ndim,Nsnap,Ncomp,Ngrid,Npod,Nfilt,f,InputSize,W,boundary,corr_type] = ...
                                           parse_input(varargin,nargin);
 
if corr_type==1 % temporal correlation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% perform snapshot SPOD                                                   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
    % calculate temporal correlation matrix %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    R=zeros(Nsnap); % size: Nsnap x Nsnap
    for i=1:Ncomp
        U = reshape(varargin{i},Ngrid,Nsnap);
        R = R + double(U)'*spdiags(W,0,Ngrid,Ngrid)*double(U);
    end
    R = R/(Nsnap*sum(W));
    
    % calculate filtered correlation matrix S %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if boundary == 1
        % convolve periodically extended matrix with filter diagonal-matrix
        ind = 1-Nfilt:Nsnap+Nfilt;
        ind = mod(ind-1, Nsnap) + 1;
        S = conv2(R(ind,ind),diag(f),'valid');
    elseif boundary == 2 %DFT case
        r = zeros(Nsnap,1);
        for i=0:Nsnap-1
            r(i+1) = sum(diag(R,i))*f(1);
        end
        % periodic boundary conditions -> circulant matrix
        r(2:end) = r(2:end)+r(end:-1:2);
        S = toeplitz(r);
    else
        % convolve zero padded matrix with filter diagonal-matrix
        S = conv2(R,diag(f),'same');
    end
    % alterative element-wise formulation is much slower
    % idx = -nf:nf;
    % for i=1:n
    %     for j=1:n
    %         ind = sub2ind([n,n],mod(i+idx,n)+1,mod(j+idx,n)+1);
    %         S(i,j) = sum(f.*R(ind));
    %     end
    % end
 
    
    % calculate eigenvalues and eigenvectors of S %%%%%%%%%%%%%%%%%%%%%%%%%
    [V,D] = eig(S);
    [lambda,idx] = sort(diag(D),'descend'); % energy of modes
    V = V(:,idx); % normalized temporal coefficients
    
    % cut of POD modes according to rank of matrix S 
    Nrank = sum(lambda>eps(lambda(1))*max(Nsnap,Ngrid*Ncomp));
    if isempty(Npod)
        Npod = Nrank;
    elseif Npod>Nrank
        Npod = Nrank;
        warning('SPOD:RankDeficit',...
    'Correlation matrix rank is less than number of requested POD modes')
    end
    
    % compute scaled temporal coefficients %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    a = bsxfun(@times,V(:,1:Npod),sqrt(Nsnap*lambda(1:Npod))');
    
    % calculate spatial modes for all components %%%%%%%%%%%%%%%%%%%%%%%%%%
    a_proj = bsxfun(@rdivide,V(:,1:Npod),sqrt(Nsnap*lambda(1:Npod))');
    OutputSize = InputSize;
    OutputSize(Ndim) = Npod;
    mode_norm = zeros(1,Npod);
    for i=1:Ncomp
        U = reshape(varargin{i},Ngrid,Nsnap);
        Ui = U*a_proj;
        mode_norm = mode_norm + sum(bsxfun(@times,Ui.^2,W),1)/sum(W);
        % reshape and put it in the output
        varargout{i} = reshape(Ui,OutputSize);
    end
    
elseif corr_type==0  % spatial correlation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% perform spatial SPOD                                                    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
    % calculate spatial correlation matrix %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % including time shifted and mixed component correlations
    Nfilt2 = length(f); % = 2*Nfilt + 1
    Ncorr = Ngrid*Nfilt2*Ncomp;
    
    R=zeros(Ngrid,Nfilt2,Ncomp,Ngrid,Nfilt2,Ncomp);
    for k=1:Ncomp % first loop over velocity component
        Uk = bsxfun(@times,reshape(varargin{k},Ngrid,Nsnap),sqrt(W));
    for l=1:Ncomp % second loop over velocity component
        Ul = bsxfun(@times,reshape(varargin{l},Ngrid,Nsnap),sqrt(W));
        for i=-Nfilt:Nfilt % first loop over filter coefficient
        for j=-Nfilt:Nfilt % second loop over filter coefficient
            indi = i+Nfilt+1; 
            indj = j+Nfilt+1; 
            ij_shift = i-j;
            f_ij = sqrt(f(indi)*f(indj));
            if boundary==1
                % correlation with periodic boundary conditions
                R(:,indi,k,:,indj,l) = f_ij*Uk*circshift(Ul,[0 ij_shift])';
            else
                % correlation with finite signal (zero padded)
                u_lim = min([Nsnap-i,Nsnap-j,Nsnap]); 
                l_lim = max([1-i,1-j,1]);
                subi = (l_lim:u_lim)+i;
                subj = (l_lim:u_lim)+j;
                R(:,indi,k,:,indj,l) = f_ij*Uk(:,subi)*Ul(:,subj)';
            end
        end
        end
    end
    end
    clear Uk Ul
    R = reshape(R,Ncorr,Ncorr)/(Nsnap*sum(W));
 
    % calculate eigenvalues and eigenvectors of R %%%%%%%%%%%%%%%%%%%%%%%%%
    [V,D] = eig(R);
    [lambda,idx] = sort(diag(D),'descend'); % energy of modes
    V = V(:,idx); % normalized spatial modes (+ time delayed modes)
    
    % cut of POD modes according to rank of matrix R if necessary
    Nrank = sum(lambda>eps(lambda(1))*max(Nsnap,Ncorr));
    if isempty(Npod)
        Npod = Nrank;
    elseif Npod>Nrank
        Npod = Nrank;
        warning('SPOD:RankDeficit',...
    'Correlation matrix rank is less than number of requested POD modes')
    end
    
    % compute scaled temporal coefficients %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    OutputSize = InputSize;
    OutputSize(Ndim) = Npod;
    a = zeros(Nsnap,Npod);
    mode_norm = zeros(1,Npod);
    Nconvfilt = Ngrid*Nfilt2;
    for j=1:Ncomp
        U = reshape(varargin{j},Ngrid,Nsnap);
        if boundary==1
            % create periodically extended time series 
            Uext = [U(:,end-Nfilt+1:end),U,U(:,1:Nfilt)];
        else
            % create zero padded time series 
            Uext = [zeros(Ngrid,Nfilt),U,zeros(Ngrid,Nfilt)];
        end
        idx = (j-1)*Nconvfilt + (1:Nconvfilt);
        for i=1:Npod
            % pick apropriate spatial mode
            f_proj = reshape(V(idx,i),Ngrid,Nfilt2);
            % apply spatial weighting
            f_proj = bsxfun(@times,f_proj,sqrt(W/sum(W))); 
            % apply temporal weighting (filter)
            f_proj = bsxfun(@times,f_proj,sqrt(f)); 
            % flip dimensions before convolution to obtain a correlation 
            f_proj = rot90(f_proj,2); 
            % calculate temporal coefficient
            a(:,i) = a(:,i) + (conv2(Uext,f_proj,'valid'))';
        end
        
        % pick central spatial mode for the output and scale it
        idx = (j-1)*Nconvfilt + Nfilt*Ngrid + (1:Ngrid);
        scale = sqrt(sum(W)/f(Nfilt+1)./W);
        scale(W==0) = 1; % avoid singularities
        Ui = bsxfun(@times,V(idx,1:Npod),scale);
        % calculate norm of spatial modes
        mode_norm = mode_norm + sum(bsxfun(@times,Ui.^2,W),1)/sum(W);
        % reshape and put it in the output
        varargout{j} = reshape(Ui,OutputSize);
    end
end
 
% add other output variables %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [Ui,Vi,Wi,...,a,lambda,mode_norm] = SPOD(...)
varargout{Ncomp+1} = a;
varargout{Ncomp+2} = lambda;
varargout{Ncomp+3} = sqrt(mode_norm);
end
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Internal Functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
% input agrgument parser 
function [Ndim,Nsnap,Ncomp,Ngrid,Npod,Nfilt,f,InputSize,...
                          W,boundary,corr_type] = parse_input(input,Ninput)
    % determine dimensions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if Ninput<1
        error('too few input arguments')
    end
 
    % determine size of input data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    idx = 1; Ndim = 1;
    while idx == 1 || (ndims(input{idx}) == Ndim && numel(input{idx})>1)
        U = input{idx};
        Ndim = ndims(U);
        if idx == 1; 
            InputSize = size(U);
            Nsnap = InputSize(Ndim);
        else
            if size(U,Ndim) ~= Nsnap
                error('number of snapshots in additional component is not consistent with first component')
            end
        end
        idx = idx + 1;
    end
    Ncomp = idx - 1;
    Ngrid = prod(InputSize(1:Ndim-1));
    
    % check additional input atguments %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % filter legth
    if Ninput<idx || isempty(input{idx})
        Nfilt = 0;
    else
        Nfilt = round(input{idx});
        if isscalar(Nfilt)
            if Nfilt>Nsnap
                warning('SPOD:NfiltToLarge',...
                    'filter size is larger than number of snapshots')
            elseif Nfilt<0
                error('filter size must be a positive number')
            end
        else
            error('filter size must scalar number')
        end
    end
    idx = idx + 1;
    
    % number of modes in output
    if Ninput<idx || isempty(input{idx})
        Npod = [];
    else
        Npod = round(input{idx});
        if isscalar(Nfilt)
            if Npod<1
                error('number of modes must be larger than zero')
            end
        else
            error('number of modes must scalar number')
        end
    end
    idx = idx + 1;
    
    % spatial weighting for inner product
    if Ninput<idx || isempty(input{idx})
        W = ones(Ngrid,1);
    else
        W = input{idx};
        if numel(W) ~= Ngrid;
           error('dimension of weighting does not match the velocity data')
        end
        W = W(:);
    end
    idx = idx + 1;
    
    % time series boundary treatment
    if Ninput<idx || isempty(input{idx})
        boundary = 1;
    else
        boundary = input{idx};
        if ~isscalar(boundary) || (boundary~=1 && boundary~=0)
            error('boundary parameter must be 0 or 1')
        end
    end
    
    % temporal or spatial correlation?
    idx = idx + 1;
    if Ninput<idx || isempty(input{idx})
        corr_type = 1;
    else
        corr_type = input{idx};
        if ~isscalar(corr_type) || (corr_type~=1 && corr_type~=0)
            error('correlation type must be 0 or 1')
        end
    end
    
    % define filter %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if Nfilt==Nsnap && corr_type==1 % only for snapshot SPOD
        % choose box filter to obtain Fourier transform 
        f = filter_coefficients(Nfilt,'box');
        boundary = 2;
    else % choose gauss filter for smooth results 
        f = filter_coefficients(Nfilt,'gauss');
    end
    
    % check problem size and display apropriate warnig %%%%%%%%%%%%%%%%%%%%
    if corr_type==1 % snapshot
        Ncorr = Nsnap;
    else % spatial
        Ncorr = Ngrid*length(f);
    end
    if Ncorr>10000
        warning('SPOD:LargeProblem','computation takes long time >1 hour')
    elseif Ncorr>2000
        warning('SPOD:LargeProblem','computation takes some time >1 minute')
    elseif Ncorr>1000
        warning('SPOD:LargeProblem','computation may take some time')
    end
    
    % check number of SPOD modes %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if ~isempty(Npod)
        if Npod>Ncorr && corr_type==1
            error('number of SPOD modes must be less than number of snapshots for snapshot SPOD')
        elseif Npod>(Ncorr*Ncomp) && corr_type==0
            error('number of SPOD modes must be less than number of gridpoints times filter size for spatial SPOD')
        end
    end
end
 
function  f = filter_coefficients(Nfilt,type)
 
switch type
    case 'box'
        f = ones(Nfilt,1);
    case 'gauss'
        f = exp(-linspace(-2.285,2.285,2*Nfilt+1).^2);
    otherwise
        error('unknown filter type')
end
f = f/sum(f(:));
 
% filter size specification:
% adjust gauss filter to same cut off frequency as the half size box filter
% Nfilt = 10;
% lim = 2.285;
% f1 = exp(-linspace(-lim,lim,2*Nfilt+1).^2);
% f2 = ones(Nfilt,1);
% f1 = f1/sum(f1(:));
% f2 = f2/sum(f2(:));
% 
% figure
% hold on
% plot(abs(fftshift(fft(f1,10000))))
% plot(abs(fftshift(fft(f2,10000))),'r')
end

