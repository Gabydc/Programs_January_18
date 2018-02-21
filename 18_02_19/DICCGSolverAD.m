classdef DICCGSolverAD < LinearSolverAD
    %Preconditioned CG solver.
    %
    % SYNOPSIS:
    %   solver = PCG_ICSolverAD()
    %
    % DESCRIPTION:
    %   Solve a linearized problem using a CG solver with IC
    %   preconditioner.
    %
    % REQUIRED PARAMETERS:
    %   None
    %
    % OPTIONAL PARAMETERS (supplied in 'key'/value pairs ('pn'/pv ...)):
    %   See class properties.
    %
    %
    % SEE ALSO:
    %   BackslashSolverAD, CPRSolverAD, LinearizedProblem, LinearSolverAD
    
    %{
Copyright 2009-2016 SINTEF ICT, Applied Mathematics.

This file is part of The MATLAB Reservoir Simulation Toolbox (MRST).

MRST is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

MRST is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with MRST.  If not, see <http://www.gnu.org/licenses/>.
    %}
    properties
        % Type: {'nofill'}, 'crout', 'ictp'
        ictype
        % Drop tolerance for elements (default: 0)
        dropTolerance
        % Modified ic type: 'row', 'col', {'off'}
        michol
        % Boolean indicating replacement of zero diagonal entries
        diagcomp
        % Pivoting threshold for ilupt. Default 1.
        shape
        %  Determines which triangle is referenced and returned
        x0
        % Initial condition
        cn
        %Compute the condition number
        W
        %If wells are included
        Z
        %Deflation matrix
        dir
        %directory to save figures
        bc
        
    end
    methods
        function solver = DICCGSolverAD(varargin)
            solver = solver@LinearSolverAD();
            %      type     --- Type of factorization.
            %      droptol  --- Drop tolerance when type is 'ict'.
            %      michol   --- Indicates whether to perform Modified incomplete Cholesky.
            %      diagcomp --- Perform Compensated incomplete Cholesky with the specified
            %                   coefficient.
            %      shape    --- Determines which triangle is referenced and returned. The default value is 'lower'
            %
            solver.ictype               = 'nofill';
            solver.dropTolerance         = 0;
            solver.michol = 'off';
            solver.diagcomp     = 0;
            solver.shape        = 'lower';
            solver.x0 = [];
            solver.W= [];
            solver.cn = 0;
            solver.Z= [];
            solver.dir=[];
            solver.bc = [];
            solver = merge_options(solver, varargin{:});
        end
        
        function [result, report] = solveLinearSystem(solver, A, b)
            
            nel = size(A, 1);
            
            [L] = ichol(A, solver.getOptsIC());
%             % Display message of options
%             dopts      = opts{1}(1);
%             % Compute condition number of the matrix A
%             A_cn       = opts{1}(2);
%             % Compute true solution
%             x_true     = opts{1}(3);
%             % Checks the convergence of the method, in the case of DICCG the residual
%             % can increas, in that case, the solution will be the solution with
%             % minimal residual
%             Convergence = opts{1}(4);
%             % Save the variables to plot residual
%             Residual    = opts{1}(5);
%             % Display the number of iterations
%             Iter_m      = opts{1}(6);
%             % Compute eigenvalues of matrix A
%             Amatrix_eigs = opts{1}(7);
%             % Compute eigenvalues of matrix M^{-1}A
%             MAmatrix_eigs = opts{1}(8);
%             % Compute eigenvalues of matrix PM^{-1}A
%             PMAmatrix_eigs = opts{1}(9);
%             % Compute eigenvalues and condition number of matrix E
%             E_cn = opts{1}(10);

            opts = {{true, false, false, true, true, true, false, false, false,false}};
            
            [result,flag,res,its,resvec,resulte]  = DICCG_MRST(A,b,solver.Z,...
                solver.tolerance,min(solver.maxIterations,nel),L,L',solver.x0, ...
                'opts',opts);
            
            report = struct('CGflag',  flag, ...
                'relres',   res,...
                'res',   resvec,...
                'iter', its, ...
                'extras',resulte);
            
            
        end
        
        function opts = getOptsIC(solver)
            opts = struct('type',    solver.ictype, ...
                'droptol', solver.dropTolerance, ...
                'michol',    solver.michol , ...
                'diagcomp',   solver.diagcomp, ...
                'shape',  solver.shape);
        end
        
        function [dx, result, report] = solveLinearProblem(solver, problem, model)
            % Reduce to cell variables before solving
            [dx, result, report]= solver.solveCellReducedLinearProblem(problem, model);
        end
    end
end
