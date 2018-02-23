function [varargout] = matfem (xe, bc_x0, bc_xL,varargin)

%==========================================================================
% The MatMol Group (2009)
%
% [fem_im] = matfem (xe, bc_x0, bc_xL, fem_om)
%
% This function is a modification of functions "sdl_sys", "udl_sys",
% "sdq_sys" and "udq_sys" included in the fselib Library (see
% http://dehesa.freeshell.org/FSELIB). This function allows us to extract
% the mass, diffusion, convection and boundary matrices using the finite
% element method. The basis functions of the FEM can be linear or
% quadratic. The method and the original functions "sdl_sys", "udl_sys",
% "sdq_sys", "udq_sys" are described in the book:
% Pozrikidis, C. (2005). Introduction to Finite and Spectral Element
% Methods using Matlab. Chapman & Hall/CRC.
%
% Input parameters:
% xe:       Coordinates of the spatial discretization points.
% bc_x0:    Type of boundary conditions at the beginning of the domain.
% bc_xL:    Type of boundary conditions at the end of the domain.
% fem_om:   This variable indicates the FEM matrices that the user wants to
%           compute. The matrix is indicated by means of a string: 'MM' for
%           the mass matrix, 'DM' for the diffusion matrix, 'CM' for the
%           convection matrix and 'BM' for the boundary matrix. The last
%           element of this variable indicates the type of basis functions
%           employed: 'lin' for linear and 'quad' for quadratic.  In the
%           case of linear basis functions the last element can be omitted.
%       
% The boundary type can be 'dir' (Dirichlet boundary conditions), 'neu'
% (Neumann boundary conditions) or 'rob' (Robin boundary conditions)
%
% Output parameters:
% fem_im:   The output parameters are the matrices indicated in fem_om.
%           They must be sorted following the same sequence used in fem_om.
%           If quadratic elements are chosen, it may result convinient to 
%           recover the mesh with the interior points. This is done by 
%           adding a new output variable.
%
% Example: Compute the mass matrix and the boundary matrix with Dirichlet
% boundary conditions in the first point and Neumann in the last point
% xe = linspace(0,pi,21);
% [mass, boundary] = matfem (xe, 'dir', 'neu', 'MM', 'BM')
%
% Example: Compute the convection matrix and the diffusion matrix with Robin
% boundary conditions in the first point and Dirichlet in the last point 
% using second order basis function
% xe = linspace(0,pi,21);
% [conv, diffus,xex] = matfem (xe, 'rob', 'dir', 'CM', 'DM','quad')
%==========================================================================

% Number of matrices requested
nip = length(varargin);
switch lower(varargin{nip})
    case {'quad','lin'}
        basis       = varargin{nip};
        n_matrices	= nip - 1;
    otherwise
        basis       = 'lin';
        n_matrices	= nip;
end


% Selection of the method
switch basis
    case{'lin'}

        %-------------
        % element size
        %-------------
        ne	= length(xe) - 1;
        for l=1:ne
            h(l) = xe(l+1)-xe(l);
        end

        %==========================================================================
        %       COMPUTATION OF THE MASS MATRIX

        %----------------------------------
        % initialize the tridiagonal matrix
        %----------------------------------
        at = zeros(ne+1,1);
        bt = zeros(ne+1,1);
        ct = zeros(ne+1,1);

        %----------------------------------
        % loop over the first ne-1 elements
        %----------------------------------
        for l=1:ne

            B11 = h(l)/3.0; B12 = 0.5*B11;  % mass matrix
            B21 = B12;      B22 = B11;

            at(l)     = at(l) + B11;        % tridiagonal matrix components
            bt(l)     = bt(l) + B12;
            ct(l)     = ct(l) + B21;
            at(l+1)   = at(l+1) + B22;

        end

        %------------------------
        % last element is special
        %------------------------
        B11 = h(ne)/3.0; B12 = 0.5*B11;

        at(ne+1) = B11;

        % Construction of the mass matrix from the diagonal vectors
        m_mat = diag(at) + diag(bt(1:end-1),1) + diag(ct(1:end-1),-1);
        clear B11 B21 B22 B12
        %==================================================================





        %==================================================================
        %       COMPUTATION OF THE DIFFUSION MATRIX

        %----------------------------------
        % initialize the tridiagonal matrix
        %----------------------------------
        at = zeros(ne+1,1);
        bt = zeros(ne+1,1);
        ct = zeros(ne+1,1);

        %----------------------------------
        % loop over the first ne elements
        %----------------------------------
        for l=1:ne

            A11 = 1/h(l); A12 =-A11;    % Diffusion matrix
            A21 = A12;    A22 = A11;

            at(l) = at(l) + A11;        % tridiagonal matrix components
            bt(l) = bt(l) + A12;
            ct(l) = ct(l) + A21;
            at(l+1) = at(l+1) + A22;

        end

        %----------------------------
        % the last element is special
        %----------------------------
        A11 = 1.0/h(ne);
        at(ne+1) = A11;

        % Construction of the diffusion matrix from the diagonal vectors
        d_mat = diag(at) + diag(bt(1:end-1),1) + diag(ct(1:end-1),-1);
        clear A11 A21 A22 A12
        %==================================================================





        %==================================================================
        %       COMPUTATION OF THE CONVECTION MATRIX

        %-----------
        % initialize
        %-----------
        at = zeros(ne+1,1);
        bt = zeros(ne+1,1);
        ct = zeros(ne+1,1);

        cf = 0.5;

        C11 = -cf; C12 = cf;        % Convection matrix
        C21 = -cf; C22 = cf;

        %------------------------
        % loop over ne elements
        %------------------------
        for l=1:ne
            at(l)   = at(l) + C11;  % tridiagonal matrix components
            bt(l)   = bt(l) + C12;
            ct(l)   = ct(l) + C21;
            at(l+1) = at(l+1) + C22;
        end

        %-------------
        % last element
        %-------------
        at(ne+1) = -C11;

        % Construction of the convection matrix from the diagonal vectors
        c_mat   = diag(at) + diag(bt(1:end-1),1) + diag(ct(1:end-1),-1);
        clear at bt ct C11 C12 C21 C22
        %==================================================================





        %==================================================================
        %       COMPUTATION OF THE BOUNDARY MATRIX
        %----------------------------------
        % initialize the boundary matrix
        %----------------------------------
        b_mat   = zeros(ne+1,ne+1);

        % Boundary conditions at the beginning of the spatial domain
        switch bc_x0
            case {'dir'}
                b_mat(1,1)  = 1000;
            case {'neu'}
                b_mat(1,1)  = 0;
            case {'rob'}
                b_mat(1,1)  = 1;
        end
        % Boundary conditions at the end of the spatial domain
        switch bc_xL
            case {'dir'}
                b_mat(ne+1, ne+1)   = 1000;
            case {'neu'}
                b_mat(ne+1, ne+1)   = 0;
            case {'rob'}
                b_mat(ne+1, ne+1)   = 1;
        end
        %==================================================================
        % OUTPUT VARIABLES
        for ii = 1 : n_matrices
            switch varargin{ii}
                case{'MM'}
                    varargout{ii} = m_mat;
                case{'DM'}
                    varargout{ii} = d_mat;
                case{'CM'}
                    varargout{ii} = c_mat;
                case{'BM'}
                    varargout{ii} = b_mat;
            end
        end





    case{'quad'}

        %-------------
        % element size
        %-------------
        ne  = length(xe) - 1;
        for l=1:ne
            h(l) = xe(l+1)-xe(l);
        end
        %------------------------------
        % number of unique global nodes
        %------------------------------
        ng  = 2*ne+1;
        for ii = 1 : ne
            xex(2*(ii-1)+1 ,1)  = xe(ii);
            xex(2*ii , 1)       = xe(ii) + (xe(ii+1)-xe(ii))/2;
        end
        xex(ng) = xe(ne+1);


        %==========================================================================
        %       COMPUTATION OF THE MASS MATRIX

        %------------------------------------
        % initialize the pentadiagonal matrix
        %------------------------------------
        ap = zeros(ng,1); bp = zeros(ng,1); cp = zeros(ng,1);
        dp = zeros(ng,1); ep = zeros(ng,1);

        %-----------------------
        % loop over all elements
        %-----------------------
        for l=1:ne

            cf = h(l)/30.0;
            B11 = 4.0*cf; B12 =  2.0*cf;  B13 = -cf;
            B21 = B12;    B22 = 16.0*cf;  B23 = B12;
            B31 = B13;    B32 = B23;      B33 = B11;

            cl1 = 2*l-1; cl2 = 2*l; cl3 = 2*l+1;

            ap(cl1) = ap(cl1) + B11;
            bp(cl1) = bp(cl1) + B12;
            cp(cl1) = cp(cl1) + B13;
            dp(cl2) = dp(cl2) + B21;
            ap(cl2) = ap(cl2) + B22;
            bp(cl2) = bp(cl2) + B23;
            ep(cl3) = ep(cl3) + B31;
            dp(cl3) = dp(cl3) + B32;
            ap(cl3) = ap(cl3) + B33;
        end

        % Construction of the diffusion matrix from the diagonal vectors
        m_mat = diag(ap) + diag(bp(1:end-1),1) + diag(cp(1:end-2),2) + diag(dp(2:end),-1) + diag(ep(3:end),-2);
        clear B11 B21 B22 B12 ap bp cp dp ep
        %==========================================================================




        %==========================================================================
        %       COMPUTATION OF THE DIFFUSION MATRIX

        %------------------------------------
        % initialize the pentadiagonal matrix
        %------------------------------------
        ap = zeros(ng,1); bp = zeros(ng,1); cp = zeros(ng,1);
        dp = zeros(ng,1); ep = zeros(ng,1);

        %-----------------------
        % loop over all elements
        %-----------------------
        for l=1:ne

            cf = 1.0/(3.0*h(l));
            A11 = 7.0*cf;  A12 = -8.0*cf; A13 = cf;
            A21 = A12;     A22 = 16.0*cf; A23 = A12;
            A31 = A13;     A32 = A23;     A33 = A11;

            cl1 = 2*l-1; cl2 = 2*l; cl3 = 2*l+1;

            ap(cl1) = ap(cl1) + A11;
            bp(cl1) = bp(cl1) + A12;
            cp(cl1) = cp(cl1) + A13;
            dp(cl2) = dp(cl2) + A21;
            ap(cl2) = ap(cl2) + A22;
            bp(cl2) = bp(cl2) + A23;
            ep(cl3) = ep(cl3) + A31;
            dp(cl3) = dp(cl3) + A32;
            ap(cl3) = ap(cl3) + A33;
        end

        % Construction of the diffusion matrix from the diagonal vectors
        d_mat = diag(ap) + diag(bp(1:end-1),1) + diag(cp(1:end-2),2) + diag(dp(2:end),-1) + diag(ep(3:end),-2);
        clear A11 A21 A22 A12 ap bp cp dp ep
        %==========================================================================




        %==========================================================================
        %       COMPUTATION OF THE DIFFUSION MATRIX

        %------------------------------------
        % initialize the pentadiagonal matrix
        %------------------------------------
        ap = zeros(ng,1); bp = zeros(ng,1); cp = zeros(ng,1);
        dp = zeros(ng,1); ep = zeros(ng,1);

        %-----------------------
        % loop over all elements
        %-----------------------
        for l=1:ne

            cf = 1.0/(6);
            C11 = -3*cf;  C12 = 4*cf; C13 = -cf;
            C21 = -C12;   C22 = 0;    C23 = C12;
            C31 = -C13;   C32 = -C23; C33 = -C11;

            cl1 = 2*l-1; cl2 = 2*l; cl3 = 2*l+1;

            ap(cl1) = ap(cl1) + C11;
            bp(cl1) = bp(cl1) + C12;
            cp(cl1) = cp(cl1) + C13;
            dp(cl2) = dp(cl2) + C21;
            ap(cl2) = ap(cl2) + C22;
            bp(cl2) = bp(cl2) + C23;
            ep(cl3) = ep(cl3) + C31;
            dp(cl3) = dp(cl3) + C32;
            ap(cl3) = ap(cl3) + C33;
        end

        % Construction of the diffusion matrix from the diagonal vectors
        c_mat = diag(ap) + diag(bp(1:end-1),1) + diag(cp(1:end-2),2) + diag(dp(2:end),-1) + diag(ep(3:end),-2);
        clear C11 C21 C22 C12 ap bp cp dp ep
        %==========================================================================

        %==================================================================
        %       COMPUTATION OF THE BOUNDARY MATRIX
        %----------------------------------
        % initialize the boundary matrix
        %----------------------------------
        b_mat   = zeros(2*ne+1,2*ne+1);

        % Boundary conditions at the beginning of the spatial domain
        switch bc_x0
            case {'dir'}
                b_mat(1,1)  = 1000;
            case {'neu'}
                b_mat(1,1)  = 0;
            case {'rob'}
                b_mat(1,1)  = 1;
        end
        % Boundary conditions at the end of the spatial domain
        switch bc_xL
            case {'dir'}
                b_mat(2*ne+1, 2*ne+1)   = 1000;
            case {'neu'}
                b_mat(2*ne+1, 2*ne+1)   = 0;
            case {'rob'}
                b_mat(2*ne+1, 2*ne+1)   = 1;
        end
        %==================================================================



        % OUTPUT VARIABLES
        for ii = 1 : n_matrices
            switch varargin{ii}
                case{'MM'}
                    varargout{ii} = m_mat;
                case{'DM'}
                    varargout{ii} = d_mat;
                case{'CM'}
                    varargout{ii} = c_mat;
                case{'BM'}
                    varargout{ii} = b_mat;
            end
        end

        % The extended mesh
        varargout{ii+1} = xex;
        if (length(varargout) > nargout)
            fprintf('Warning: It may be convenient to add a new output variable in order to recover the extended mesh \n')
        end
end


return;
