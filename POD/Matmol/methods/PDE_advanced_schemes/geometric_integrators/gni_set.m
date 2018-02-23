%...  The MatMol Group (2016)
    function options = gni_set(varargin)

% GNISET Create/alter GNI OPTIONS structure.
%   The syntax for GNI_SET is the same as for ODESET, but the list of possible
%   properties is different. 
%   OPTIONS = GNI_SET('NAME1',VALUE1,'NAME2',VALUE2,...) creates an integrator
%   options structure OPTIONS in which the named properties have the
%   specified values. Any unspecified properties have default values.
%   
%   OPTIONS = GNI_SET(OLDOPTS,'NAME1',VALUE1,...) alters an existing options
%   structure OLDOPTS.
%   
%   OPTIONS = GNI_SET(OLDOPTS,NEWOPTS) combines an existing options structure
%   OLDOPTS with a new options structure NEWOPTS. Any new properties
%   overwrite corresponding old properties. 
%
% GNISET PROPERTIES
% StepSize - Step size  [ positive scalar ]
%   For solver that will use a fixed step size. The step size is given by 
%   StepSize. It may be slightly altered if the length of the integration
%   interval is not an integer multiple of the step size.
%  
% NumSteps - Number of steps  [ positive integer ]
%   The solver will use a fixed step size adjusted to make NumSteps steps.
%   If StepSize is set, it overrides this option.
%   
% Precision - precision [positive scalar]
%   Parameter of precision for variable step size solver
%
% MaxIter - Maximal number of iterations [ integer | {50} ]
%   Maximal number of fixed point iterations performed at each step for 
%   gni_fgauss1.
%
% Method - Name of the method [ string | {'G12'} ]
%   For gni_fgauss1 and gni_vgauss1, implemented methods are 'G2', 'G4', 
%       'G6', 'G8', 'G10', 'G12' and 'G14'. 
%   For gni_fstrvl2 and gni_vstrvl2, implemented methods are '21', '43',
%       '45', '67', '69', '815', '817', '1031', '1033', and '1035'. 
%   for gni_flmm2, implemented methods are '4a', '4b', '6a', '6b', '8a',
%       '8b', '8c', '10a', '10b', and '12a'. 
%   In all cases the first one or two digits denote the order of the
%   method.
%
% Vectorized - Vectorized ODE file  [ on | {off} ]
%   Set this property 'on' if the ODE file is coded so that F(t,[y1 y2 ...])
%   returns [F(t,y1) F(t,y2) ...].
%   
% RescalFun - Name of the rescaling function [ string ]
%   The rescaling function is a positive scalar function used by the
%   stormer-verlet method gni_vstrvl2.
%
% OutputFcn - Installable output function  [ string ]
%   This output function is called by the solver after each time step. When
%   a solver is called with no output arguments, OutputFcn defaults to 
%   'odeplot'. Otherwise, OutputFcn defaults to [].
%
% OutputSel - Output selection indices  [ vector of integers ]
%   This vector of indices specifies which components of the solution vector
%   are passed to the OutputFcn. OutputSel defaults to all components.
%   
% OutputSteps - Which steps to output  [ integer | 1 ]
%   This value tells which computed solution points are sent to the output. If
%   OutputSteps = 10, every 10th solution point is sent to the output. 
%  
% Print out possible values of properties.
if (nargin == 0) && (nargout == 0)
  fprintf('        StepSize: [ positive scalar    ]\n');
  fprintf('        NumSteps: [ positive integer   ]\n');
  fprintf('       Precision: [ positive scalar    ]\n');
  fprintf('         MaxIter: [ positive integer   ]\n');
  fprintf('          Method: [ string             ]\n');
  fprintf('      Vectorized: [ on | {off}         ]\n');
  fprintf('       RescalFun: [ string             ]\n');
  fprintf('     Compensated: [ on | {off}         ]\n');
  fprintf('        Momentum: [ on | {off}         ]\n');
  fprintf('       OutputFcn: [ string             ]\n');
  fprintf('       OutputSel: [ vector of integers ]\n');
  fprintf('     OutputSteps: [ positive integer   ]\n');
  fprintf('\n');
  return;
end

Names = {
	'StepSize'
	'NumSteps'
    'Precision'
	'MaxIter'
	'Method'
    'Vectorized'
    'RescalFun'
    'Compensated'
    'Momentum'
    'OutputFcn'
    'OutputSel'
    'OutputSteps'
    };
m = size(Names,1);

options = [];

i = 1;

while i<=nargin
    arg = varargin{i};
    if ischar(arg)
        break
    end   
    if ~isempty(arg)
        if ~isa(arg,'struct')
            error(['Expected argument %d to be a string property name ' ...
                'or an options structure \n created with GNI_SET.'],i);
        end        
        if isempty(options)
            options = arg;
        else
            for j = 1:m
                val = arg.(strjoin(Names(j,:)));
                if ~isequal(val,[])  
                    options.(strjoin(Names(j,:))) = val;                   
                end
            end
        end
    end    
    i = i + 1;
end

if isempty(options)
    for j = 1:m
        options.(strjoin(Names(j,:))) = [];
    end
end

if rem(nargin-i+1,2) ~= 0
    error('Arguments must occur in name-value pairs.');
end

expectval = 0;

while i <= nargin
    arg = varargin{i};    
    if ~expectval
        if ~ischar(arg)
            error('Expected argument %d to be a string property name.',i);
        end
        lowArg = lower(arg);
        j = strncmpi(Names,lowArg,length(lowArg));
        l = find(j);
        
        if max(j) == 0
            error('Unrecognized property name ''%s''.', arg);
        elseif length(l) > 1
            k = strcmpi(lowArg,Names);  
            if length(k) == 1
                l = find(k);
            else
                msg = sprintf('Ambiguous property name ''%s'' ', arg);
                msg = [msg '(' strjoin(Names(l(1),:))];
                for k = 2:length(l)
                    msg = [msg ', ' strjoin(Names(l(k),:))];
                end
                msg = sprintf('%s).', msg);
                error(msg);
            end
        end      
        expectval = 1;
    else 
        options.(strjoin(Names(l,:))) = arg;
        expectval = 0;
    end
    i = i + 1;
end