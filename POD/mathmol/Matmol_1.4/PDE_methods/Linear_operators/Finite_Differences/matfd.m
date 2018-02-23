function [Der] = matfd (varargin)

%==========================================================================
% The MatMol Group (2009)
%
% [Deriv] = matfd (z, grid, der_ord, ns, mth, v)
%
% This function allows us to extract the finite differences matrices
%
% Input parameters:
% z:        Coordinates of the spatial discretization points.
% grid:     Set this parameter to 'uni' if the grid is uniform, otherwise,
%           set it to 'non_uni'.
% der_ord:  Order of the derivative. The possibilities are '1st', '2nd',
%           '3rd'.
% ns:       Number of points in the stencil.
% mth:      Finite differences scheme. The possibilities are 'centered',
%           'upwind', 'upwind_biased'.
% v:        Fluid velocity. It is necessary in the cases of upwind and
%           upwind biased methods
%
% Output parameters:
% Deriv:    Finite difference approximation of the derivative operator
%           of order indicated in the input parameter der_ord.
%
% Note:     Some of the different possibilities were not included either
%           because they do not exist (like the case of a centered scheme
%           with four points in the stencil) or because the authors decided
%           not to include them in this version of the MATMOL toolbox. 
%           Future versions will include new possibilities. For checking 
%           all the available possibilities, the user is refferred to the 
%           manual of this function.
%
% Example: Compute the first order finite difference matrix with uniform
%          grid using 7 points in the stencil and a centered finite 
%          difference scheme
%
%          xe   = linspace(0,pi,21);
%          [D1] = matfd (xe, 'uni', '1st', 7, 'centered')
%==========================================================================

% Input parameters
z       = varargin{1};
grid    = varargin{2};
der_ord = varargin{3};
ns      = varargin{4};
mth     = varargin{5};

% Check the number of points in the stencil
switch ns
    case {11}
        % Check the uniformity of the grid
        switch grid
            case {'uni'}
                % Check the order of the derivative
                switch der_ord
                    case {'1st'}
                        % Check the method
                        switch mth
                            case {'centered'}
                                z0    = z(1);  zL = z(end); n = length(z);
                                [Der] = eleven_point_centered_uni_D1(z0,zL,n);
                            case {'upwind'}
                                fprintf('\n This combination is not possible in the version of the MATMOL.\n')
                                return
                            case {'upwind_biased'}
                                fprintf('\n This combination is not possible in the version of the MATMOL.\n')
                                return
                            otherwise
                                fprintf('\n Unknown method.\n')
                        end
                    case {'2nd'}
                        % Check the method
                        switch mth
                            case {'centered'}
                                fprintf('\n This combination is not possible in the version of the MATMOL.\n')
                                return
                            case {'upwind'}
                                fprintf('\n This combination is not possible in the version of the MATMOL.\n')
                                return
                            case {'upwind_biased'}
                                fprintf('\n This combination is not possible in the version of the MATMOL.\n')
                                return
                            otherwise
                                fprintf('\n Unknown method.\n')
                        end
                    case {'3rd'}
                        % Check the method
                        switch mth
                            case {'centered'}
                                fprintf('\n This combination is not possible in the version of the MATMOL.\n')
                                return
                            case {'upwind'}
                                fprintf('\n This combination is not possible in the version of the MATMOL.\n')
                                return
                            case {'upwind_biased'}
                                fprintf('\n This combination is not possible in the version of the MATMOL.\n')
                                return
                            otherwise
                                fprintf('\n Unknown method.\n')
                        end
                    otherwise
                        fprintf('\n The maximum order of the derivative is 3.\n')
                end
            case {'non_uni'}
                % Check the order of the derivative
                switch der_ord
                    case {'1st'}
                        % Check the method
                        switch mth
                            case {'centered'}
                                [Der] = eleven_point_centered_D1(z);
                            case {'upwind'}
                                fprintf('\n This combination is not possible in the version of the MATMOL.\n')
                                return
                            case {'upwind_biased'}
                                fprintf('\n This combination is not possible in the version of the MATMOL.\n')
                                return
                            otherwise
                                fprintf('\n Unknown method.\n')
                        end
                    case {'2nd'}
                        % Check the method
                        switch mth
                            case {'centered'}
                                fprintf('\n This combination is not possible in the version of the MATMOL.\n')
                                return
                            case {'upwind'}
                                fprintf('\n This combination is not possible in the version of the MATMOL.\n')
                                return
                            case {'upwind_biased'}
                                fprintf('\n This combination is not possible in the version of the MATMOL.\n')
                                return
                            otherwise
                                fprintf('\n Unknown method.\n')
                        end
                    case {'3rd'}
                        % Check the method
                        switch mth
                            case {'centered'}
                                fprintf('\n This combination is not possible in the version of the MATMOL.\n')
                                return
                            case {'upwind'}
                                fprintf('\n This combination is not possible in the version of the MATMOL.\n')
                                return
                            case {'upwind_biased'}
                                fprintf('\n This combination is not possible in the version of the MATMOL.\n')
                                return
                            otherwise
                                fprintf('\n Unknown method.\n')
                        end
                    otherwise
                        fprintf('\n The maximum order of the derivative is 3.\n')
                end
            otherwise
                fprintf('\n The grid option is not valid.\n')
        end

    case {9}
        % Check the uniformity of the grid
        switch grid
            case {'uni'}
                % Check the order of the derivative
                switch der_ord
                    case {'1st'}
                        % Check the method
                        switch mth
                            case {'centered'}
                                z0    = z(1);  zL = z(end); n = length(z);
                                [Der] = nine_point_centered_uni_D1(z0,zL,n);
                            case {'upwind'}
                                fprintf('\n This combination is not possible in the version of the MATMOL.\n')
                                return
                            case {'upwind_biased'}
                                fprintf('\n This combination is not possible in the version of the MATMOL.\n')
                                return
                            otherwise
                                fprintf('\n Unknown method.\n')
                        end
                    case {'2nd'}
                        % Check the method
                        switch mth
                            case {'centered'}
                                fprintf('\n This combination is not possible in the version of the MATMOL.\n')
                                return
                            case {'upwind'}
                                fprintf('\n This combination is not possible in the version of the MATMOL.\n')
                                return
                            case {'upwind_biased'}
                                fprintf('\n This combination is not possible in the version of the MATMOL.\n')
                                return
                            otherwise
                                fprintf('\n Unknown method.\n')
                        end
                    case {'3rd'}
                        % Check the method
                        switch mth
                            case {'centered'}
                                fprintf('\n This combination is not possible in the version of the MATMOL.\n')
                                return
                            case {'upwind'}
                                fprintf('\n This combination is not possible in the version of the MATMOL.\n')
                                return
                            case {'upwind_biased'}
                                fprintf('\n This combination is not possible in the version of the MATMOL.\n')
                                return
                            otherwise
                                fprintf('\n Unknown method.\n')
                        end
                    otherwise
                        fprintf('\n The maximum order of the derivative is 3.\n')
                end
            case {'non_uni'}
                % Check the order of the derivative
                switch der_ord
                    case {'1st'}
                        % Check the method
                        switch mth
                            case {'centered'}
                                [Der] = nine_point_centered_D1(z);
                            case {'upwind'}
                                fprintf('\n This combination is not possible in the version of the MATMOL.\n')
                                return
                            case {'upwind_biased'}
                                fprintf('\n This combination is not possible in the version of the MATMOL.\n')
                                return
                            otherwise
                                fprintf('\n Unknown method.\n')
                        end
                    case {'2nd'}
                        % Check the method
                        switch mth
                            case {'centered'}
                                fprintf('\n This combination is not possible in the version of the MATMOL.\n')
                                return
                            case {'upwind'}
                                fprintf('\n This combination is not possible in the version of the MATMOL.\n')
                                return
                            case {'upwind_biased'}
                                fprintf('\n This combination is not possible in the version of the MATMOL.\n')
                                return
                            otherwise
                                fprintf('\n Unknown method.\n')
                        end
                    case {'3rd'}
                        % Check the method
                        switch mth
                            case {'centered'}
                                fprintf('\n This combination is not possible in the version of the MATMOL.\n')
                                return
                            case {'upwind'}
                                fprintf('\n This combination is not possible in the version of the MATMOL.\n')
                                return
                            case {'upwind_biased'}
                                fprintf('\n This combination is not possible in the version of the MATMOL.\n')
                                return
                            otherwise
                                fprintf('\n Unknown method.\n')
                        end
                    otherwise
                        fprintf('\n The maximum order of the derivative is 3.\n')
                end
            otherwise
                fprintf('\n The grid option is not valid.\n')
        end
    case {7}
        % Check the uniformity of the grid
        switch grid
            case {'uni'}
                % Check the order of the derivative
                switch der_ord
                    case {'1st'}
                        % Check the method
                        switch mth
                            case {'centered'}
                                z0    = z(1);  zL = z(end); n = length(z);
                                [Der] = seven_point_centered_uni_D1(z0,zL,n);
                            case {'upwind'}
                                fprintf('\n This combination is not possible in the version of the MATMOL.\n')
                                return
                            case {'upwind_biased'}
                                fprintf('\n This combination is not possible in the version of the MATMOL.\n')
                                return
                            otherwise
                                fprintf('\n Unknown method.\n')
                        end
                    case {'2nd'}
                        % Check the method
                        switch mth
                            case {'centered'}
                                fprintf('\n This combination is not possible in the version of the MATMOL.\n')
                                return
                            case {'upwind'}
                                fprintf('\n This combination is not possible in the version of the MATMOL.\n')
                                return
                            case {'upwind_biased'}
                                fprintf('\n This combination is not possible in the version of the MATMOL.\n')
                                return
                            otherwise
                                fprintf('\n Unknown method.\n')
                        end
                    case {'3rd'}
                        % Check the method
                        switch mth
                            case {'centered'}
                                fprintf('\n This combination is not possible in the version of the MATMOL.\n')
                                return
                            case {'upwind'}
                                fprintf('\n This combination is not possible in the version of the MATMOL.\n')
                                return
                            case {'upwind_biased'}
                                fprintf('\n This combination is not possible in the version of the MATMOL.\n')
                                return
                            otherwise
                                fprintf('\n Unknown method.\n')
                        end
                    otherwise
                        fprintf('\n The maximum order of the derivative is 3.\n')
                end
            case {'non_uni'}
                % Check the order of the derivative
                switch der_ord
                    case {'1st'}
                        % Check the method
                        switch mth
                            case {'centered'}
                                [Der] = seven_point_centered_D1(z);
                            case {'upwind'}
                                fprintf('\n This combination is not possible in the version of the MATMOL.\n')
                                return
                            case {'upwind_biased'}
                                fprintf('\n This combination is not possible in the version of the MATMOL.\n')
                                return
                            otherwise
                                fprintf('\n Unknown method.\n')
                        end
                    case {'2nd'}
                        % Check the method
                        switch mth
                            case {'centered'}
                                fprintf('\n This combination is not possible in the version of the MATMOL.\n')
                                return
                            case {'upwind'}
                                fprintf('\n This combination is not possible in the version of the MATMOL.\n')
                                return
                            case {'upwind_biased'}
                                fprintf('\n This combination is not possible in the version of the MATMOL.\n')
                                return
                            otherwise
                                fprintf('\n Unknown method.\n')
                        end
                    case {'3rd'}
                        % Check the method
                        switch mth
                            case {'centered'}
                                [Der] = seven_point_centered_D3(z);
                            case {'upwind'}
                                fprintf('\n This combination is not possible in the version of the MATMOL.\n')
                                return
                            case {'upwind_biased'}
                                fprintf('\n This combination is not possible in the version of the MATMOL.\n')
                                return
                            otherwise
                                fprintf('\n Unknown method.\n')
                        end
                    otherwise
                        fprintf('\n The maximum order of the derivative is 3.\n')
                end
            otherwise
                fprintf('\n The grid option is not valid.\n')
        end
    case {6}
        % Check the uniformity of the grid
        switch grid
            case {'uni'}
                % Check the order of the derivative
                switch der_ord
                    case {'1st'}
                        % Check the method
                        switch mth
                            case {'centered'}
                            case {'upwind'}
                            case {'upwind_biased'}
                            otherwise
                                fprintf('\n Unknown method.\n')
                        end
                    case {'2nd'}
                        % Check the method
                        switch mth
                            case {'centered'}
                            case {'upwind'}
                            case {'upwind_biased'}
                            otherwise
                                fprintf('\n Unknown method.\n')
                        end
                    case {'3rd'}
                        % Check the method
                        switch mth
                            case {'centered'}
                            case {'upwind'}
                            case {'upwind_biased'}
                            otherwise
                                fprintf('\n Unknown method.\n')
                        end
                    otherwise
                        fprintf('\n The maximum order of the derivative is 3.\n')
                end
            case {'non_uni'}
                % Check the order of the derivative
                switch der_ord
                    case {'1st'}
                        % Check the method
                        switch mth
                            case {'centered'}
                                fprintf('\n This combination is not possible in the version of the MATMOL.\n')
                                return
                            case {'upwind'}
                                fprintf('\n This combination is not possible in the version of the MATMOL.\n')
                                return
                            case {'upwind_biased'}
                                fprintf('\n This combination is not possible in the version of the MATMOL.\n')
                                return
                            otherwise
                                fprintf('\n Unknown method.\n')
                        end
                    case {'2nd'}
                        % Check the method
                        switch mth
                            case {'centered'}
                                fprintf('\n This combination is not possible in the version of the MATMOL.\n')
                                return
                            case {'upwind'}
                                fprintf('\n This combination is not possible in the version of the MATMOL.\n')
                                return
                            case {'upwind_biased'}
                                fprintf('\n This combination is not possible in the version of the MATMOL.\n')
                                return
                            otherwise
                                fprintf('\n Unknown method.\n')
                        end
                    case {'3rd'}
                        % Check the method
                        switch mth
                            case {'centered'}
                                fprintf('\n This combination is not possible in the version of the MATMOL.\n')
                                return
                            case {'upwind'}
                                fprintf('\n This combination is not possible in the version of the MATMOL.\n')
                                return
                            case {'upwind_biased'}
                                fprintf('\n This combination is not possible in the version of the MATMOL.\n')
                                return
                            otherwise
                                fprintf('\n Unknown method.\n')
                        end
                    otherwise
                        fprintf('\n The maximum order of the derivative is 3.\n')
                end
            otherwise
                fprintf('\n The grid option is not valid.\n')
        end
    case {5}
        % Check the uniformity of the grid
        switch grid
            case {'uni'}
                % Check the order of the derivative
                switch der_ord
                    case {'1st'}
                        % Check the method
                        switch mth
                            case {'centered'}
                                z0    = z(1);  zL = z(end); n = length(z);
                                [Der] = five_point_centered_uni_D1(z0,zL,n);
                            case {'upwind'}
                                fprintf('\n This combination is not possible in the version of the MATMOL.\n')
                                return
                            case {'upwind_biased'}
                                if (length(varargin) <= 5)
                                    fprintf('\n The velocity parameter is missed.\n')
                                    return
                                else
                                    v = varargin{6};
                                end
                                z0    = z(1);  zL = z(end); n = length(z);
                                [Der] = five_point_biased_upwind_uni_D1(z0,zL,n,v);
                            otherwise
                                fprintf('\n Unknown method.\n')
                        end
                    case {'2nd'}
                        % Check the method
                        switch mth
                            case {'centered'}
                                z0    = z(1);  zL = z(end); n = length(z);
                                [Der] = five_point_centered_uni_D2(z0,zL,n);
                            case {'upwind'}
                                fprintf('\n This combination is not possible in the version of the MATMOL.\n')
                                return
                            case {'upwind_biased'}
                                fprintf('\n This combination is not possible in the version of the MATMOL.\n')
                                return
                            otherwise
                                fprintf('\n Unknown method.\n')
                        end
                    case {'3rd'}
                        % Check the method
                        switch mth
                            case {'centered'}
                                fprintf('\n This combination is not possible in the version of the MATMOL.\n')
                                return
                            case {'upwind'}
                                fprintf('\n This combination is not possible in the version of the MATMOL.\n')
                                return
                            case {'upwind_biased'}
                                fprintf('\n This combination is not possible in the version of the MATMOL.\n')
                                return
                            otherwise
                                fprintf('\n Unknown method.\n')
                        end
                    otherwise
                        fprintf('\n The maximum order of the derivative is 3.\n')
                end
            case {'non_uni'}
                % Check the order of the derivative
                switch der_ord
                    case {'1st'}
                        % Check the method
                        switch mth
                            case {'centered'}
                                [Der] = five_point_centered_D1(z);
                            case {'upwind'}
                                fprintf('\n This combination is not possible in the version of the MATMOL.\n')
                                return
                            case {'upwind_biased'}
                                if (length(varargin) <= 5)
                                    fprintf('\n The velocity parameter is missed.\n'), return
                                else
                                    v     = varargin{6};
                                    [Der] = five_point_biased_upwind_D1(z,v);
                                end
                                
                            otherwise
                                fprintf('\n Unknown method.\n')
                        end
                    case {'2nd'}
                        % Check the method
                        switch mth
                            case {'centered'}
                                [Der] = five_point_centered_D2(z);
                            case {'upwind'}
                                fprintf('\n This combination is not possible in the version of the MATMOL.\n')
                                return
                            case {'upwind_biased'}
                                fprintf('\n This combination is not possible in the version of the MATMOL.\n')
                                return
                            otherwise
                                fprintf('\n Unknown method.\n')
                        end
                    case {'3rd'}
                        % Check the method
                        switch mth
                            case {'centered'}
                                fprintf('\n This combination is not possible in the version of the MATMOL.\n')
                                return
                            case {'upwind'}
                                fprintf('\n This combination is not possible in the version of the MATMOL.\n')
                                return
                            case {'upwind_biased'}
                                fprintf('\n This combination is not possible in the version of the MATMOL.\n')
                                return
                            otherwise
                                fprintf('\n Unknown method.\n')
                        end
                    otherwise
                        fprintf('\n The maximum order of the derivative is 3.\n')
                end
            otherwise
                fprintf('\n The grid option is not valid.\n')
        end
    case {4}
        % Check the uniformity of the grid
        switch grid
            case {'uni'}
                % Check the order of the derivative
                switch der_ord
                    case {'1st'}
                        % Check the method
                        switch mth
                            case {'centered'}
                                fprintf('\n This combination is not possible in the version of the MATMOL.\n')
                                return
                            case {'upwind'}
                                if (length(varargin) <= 5)
                                    fprintf('\n The velocity parameter is missed.\n'), return
                                else
                                    v = varargin{6};
                                end
                                z0    = z(1);  zL = z(end); n = length(z);
                                [Der] = four_point_upwind_uni_D1(z0,zL,n,v);
                            case {'upwind_biased'}
                                if (length(varargin) <= 5)
                                    fprintf('\n The velocity parameter is missed.\n'), return
                                else
                                    v = varargin{6};
                                    z0    = z(1);  zL = z(end); n = length(z);
                                    [Der] = four_point_biased_upwind_uni_D1(z0,zL,n,v);
                                end
                                
                            otherwise
                                fprintf('\n Unknown method.\n')
                        end
                    case {'2nd'}
                        % Check the method
                        switch mth
                            case {'centered'}
                                fprintf('\n This combination is not possible in the version of the MATMOL.\n')
                                return
                            case {'upwind'}
                                fprintf('\n This combination is not possible in the version of the MATMOL.\n')
                                return
                            case {'upwind_biased'}
                                fprintf('\n This combination is not possible in the version of the MATMOL.\n')
                                return
                            otherwise
                                fprintf('\n Unknown method.\n')
                        end
                    case {'3rd'}
                        % Check the method
                        switch mth
                            case {'centered'}
                                fprintf('\n This combination is not possible in the version of the MATMOL.\n')
                                return
                            case {'upwind'}
                                fprintf('\n This combination is not possible in the version of the MATMOL.\n')
                                return
                            case {'upwind_biased'}
                                fprintf('\n This combination is not possible in the version of the MATMOL.\n')
                                return
                            otherwise
                                fprintf('\n Unknown method.\n')
                        end
                    otherwise
                        fprintf('\n The maximum order of the derivative is 3.\n')
                end
            case {'non_uni'}
                % Check the order of the derivative
                switch der_ord
                    case {'1st'}
                        % Check the method
                        switch mth
                            case {'centered'}
                                fprintf('\n This combination is not possible in the version of the MATMOL.\n')
                                return
                            case {'upwind'}
                                if (length(varargin) <= 5)
                                    fprintf('\n The velocity parameter is missed.\n'), return
                                else
                                    v = varargin{6};
                                end
                                [Der] = four_point_upwind_D1(z,v);
                            case {'upwind_biased'}
                                if (length(varargin) <= 5)
                                    fprintf('\n The velocity parameter is missed.\n'), return
                                else
                                    v = varargin{6};
                                end
                                [Der] = four_point_biased_upwind_D1(z,v);
                            otherwise
                                fprintf('\n Unknown method.\n')
                        end
                    case {'2nd'}
                        % Check the method
                        switch mth
                            case {'centered'}
                                fprintf('\n This combination is not possible in the version of the MATMOL.\n')
                                return
                            case {'upwind'}
                                fprintf('\n This combination is not possible in the version of the MATMOL.\n')
                                return
                            case {'upwind_biased'}
                                fprintf('\n This combination is not possible in the version of the MATMOL.\n')
                                return
                            otherwise
                                fprintf('\n Unknown method.\n')
                        end
                    case {'3rd'}
                        % Check the method
                        switch mth
                            case {'centered'}
                                fprintf('\n This combination is not possible in the version of the MATMOL.\n')
                                return
                            case {'upwind'}
                                fprintf('\n This combination is not possible in the version of the MATMOL.\n')
                                return
                            case {'upwind_biased'}
                                fprintf('\n This combination is not possible in the version of the MATMOL.\n')
                                return
                            otherwise
                                fprintf('\n Unknown method.\n')
                        end
                    otherwise
                        fprintf('\n The maximum order of the derivative is 3.\n')
                end
            otherwise
                fprintf('\n The grid option is not valid.\n')
        end
    case {3}
        % Check the uniformity of the grid
        switch grid
            case {'uni'}
                % Check the order of the derivative
                switch der_ord
                    case {'1st'}
                        % Check the method
                        switch mth
                            case {'centered'}
                                z0    = z(1);  zL = z(end); n = length(z);
                                [Der] = three_point_centered_uni_D1(z0,zL,n);
                            case {'upwind'}
                                if (length(varargin) <= 5)
                                    fprintf('\n The velocity parameter is missed.\n'), return
                                else
                                    v = varargin{6};
                                end
                                z0    = z(1);  zL = z(end); n = length(z);
                                [Der] = three_point_upwind_uni_D1(z0,zL,n,v);
                            case {'upwind_biased'}
                                fprintf('\n This combination is not possible in the version of the MATMOL.\n')
                                return
                            otherwise
                                fprintf('\n Unknown method.\n')
                        end
                    case {'2nd'}
                        % Check the method
                        switch mth
                            case {'centered'}
                                z0    = z(1);  zL = z(end); n = length(z);
                                [Der] = three_point_centered_uni_D2(z0,zL,n);
                            case {'upwind'}
                                fprintf('\n This combination is not possible in the version of the MATMOL.\n')
                                return
                            case {'upwind_biased'}
                                fprintf('\n This combination is not possible in the version of the MATMOL.\n')
                                return
                            otherwise
                                fprintf('\n Unknown method.\n')
                        end
                    case {'3rd'}
                        % Check the method
                        switch mth
                            case {'centered'}
                                fprintf('\n This combination is not possible in the version of the MATMOL.\n')
                                return
                            case {'upwind'}
                                fprintf('\n This combination is not possible in the version of the MATMOL.\n')
                                return
                            case {'upwind_biased'}
                                fprintf('\n This combination is not possible in the version of the MATMOL.\n')
                                return
                            otherwise
                                fprintf('\n Unknown method.\n')
                        end
                    otherwise
                        fprintf('\n The maximum order of the derivative is 3.\n')
                end
            case {'non_uni'}
                % Check the order of the derivative
                switch der_ord
                    case {'1st'}
                        % Check the method
                        switch mth
                            case {'centered'}
                               [Der] = three_point_centered_D1(z);
                            case {'upwind'}
                                if (length(varargin) <= 5)
                                    fprintf('\n The velocity parameter is missed.\n'), return
                                else
                                    v = varargin{6};
                                end
                                [Der] = three_point_upwind_D1(z,v);
                            case {'upwind_biased'}
                                fprintf('\n This combination is not possible in the version of the MATMOL.\n')
                                return
                            otherwise
                                fprintf('\n Unknown method.\n')
                        end
                    case {'2nd'}
                        % Check the method
                        switch mth
                            case {'centered'}
                                [Der] = three_point_centered_D2(z);
                            case {'upwind'}
                                fprintf('\n This combination is not possible in the version of the MATMOL.\n')
                                return
                            case {'upwind_biased'}
                                fprintf('\n This combination is not possible in the version of the MATMOL.\n')
                                return
                            otherwise
                                fprintf('\n Unknown method.\n')
                        end
                    case {'3rd'}
                        % Check the method
                        switch mth
                            case {'centered'}
                                fprintf('\n This combination is not possible in the version of the MATMOL.\n')
                                return
                            case {'upwind'}
                                fprintf('\n This combination is not possible in the version of the MATMOL.\n')
                                return
                            case {'upwind_biased'}
                                fprintf('\n This combination is not possible in the version of the MATMOL.\n')
                                return
                            otherwise
                                fprintf('\n Unknown method.\n')
                        end
                    otherwise
                        fprintf('\n The maximum order of the derivative is 3.\n')
                end
            otherwise
                fprintf('\n The grid option is not valid.\n')
        end
    case {2}
        % Check the uniformity of the grid
        switch grid
            case {'uni'}
                % Check the order of the derivative
                switch der_ord
                    case {'1st'}
                        % Check the method
                        switch mth
                            case {'centered'}
                                fprintf('\n This combination is not possible in the version of the MATMOL.\n')
                                return
                            case {'upwind'}
                                if (length(varargin) <= 5)
                                    fprintf('\n The velocity parameter is missed.\n'), return
                                else
                                    v = varargin{6};
                                end
                                z0    = z(1);  zL = z(end); n = length(z);
                                [Der] = two_point_upwind_uni_D1(z0,zL,n,v);
                            case {'upwind_biased'}
                                fprintf('\n This combination is not possible in the version of the MATMOL.\n')
                                return
                            otherwise
                                fprintf('\n Unknown method.\n')
                        end
                    case {'2nd'}
                        % Check the method
                        switch mth
                            case {'centered'}
                                fprintf('\n This combination is not possible in the version of the MATMOL.\n')
                                return
                            case {'upwind'}
                                fprintf('\n This combination is not possible in the version of the MATMOL.\n')
                                return
                            case {'upwind_biased'}
                                fprintf('\n This combination is not possible in the version of the MATMOL.\n')
                                return
                            otherwise
                                fprintf('\n Unknown method.\n')
                        end
                    case {'3rd'}
                        % Check the method
                        switch mth
                            case {'centered'}
                                fprintf('\n This combination is not possible in the version of the MATMOL.\n')
                                return
                            case {'upwind'}
                                fprintf('\n This combination is not possible in the version of the MATMOL.\n')
                                return
                            case {'upwind_biased'}
                                fprintf('\n This combination is not possible in the version of the MATMOL.\n')
                                return
                            otherwise
                                fprintf('\n Unknown method.\n')
                        end
                    otherwise
                        fprintf('\n The maximum order of the derivative is 3.\n')
                end
            case {'non_uni'}
                % Check the order of the derivative
                switch der_ord
                    case {'1st'}
                        % Check the method
                        switch mth
                            case {'centered'}
                                fprintf('\n This combination is not possible in the version of the MATMOL.\n')
                                return
                            case {'upwind'}
                                if (length(varargin) <= 5)
                                    fprintf('\n The velocity parameter is missed.\n'), return
                                else
                                    v = varargin{6};
                                end
                                [Der] = two_point_upwind_D1(z,v);
                            case {'upwind_biased'}
                                fprintf('\n This combination is not possible in the version of the MATMOL.\n')
                                return
                            otherwise
                                fprintf('\n Unknown method.\n')
                        end
                    case {'2nd'}
                        % Check the method
                        switch mth
                            case {'centered'}
                                fprintf('\n This combination is not possible in the version of the MATMOL.\n')
                                return
                            case {'upwind'}
                                fprintf('\n This combination is not possible in the version of the MATMOL.\n')
                                return
                            case {'upwind_biased'}
                                fprintf('\n This combination is not possible in the version of the MATMOL.\n')
                                return
                            otherwise
                                fprintf('\n Unknown method.\n')
                        end
                    case {'3rd'}
                        % Check the method
                        switch mth
                            case {'centered'}
                                fprintf('\n This combination is not possible in the version of the MATMOL.\n')
                                return
                            case {'upwind'}
                                fprintf('\n This combination is not possible in the version of the MATMOL.\n')
                                return
                            case {'upwind_biased'}
                                fprintf('\n This combination is not possible in the version of the MATMOL.\n')
                                return
                            otherwise
                                fprintf('\n Unknown method.\n')
                        end
                    otherwise
                        fprintf('\n The maximum order of the derivative is 3.\n')
                end
            otherwise
                fprintf('\n The grid option is not valid.\n')
        end
    otherwise
        fprintf('\n Error: The number of points in the stencil is not allowed to be %d.\n Please see the reference manual to check the different possiblities.\n',ns)
end




return;
