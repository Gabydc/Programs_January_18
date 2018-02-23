%...  The MatMol Group (2016)
    function out = gni_get(options,name,default)

%GNI_get Get ODE OPTIONS parameters.
%   VAL = GNI_GET(OPTIONS,'NAME') extracts the value of the named property
%   from integrator options structure OPTIONS, returning an empty matrix if
%   the property value is not specified in OPTIONS. It is sufficient to type
%   only the leading characters that uniquely identify the property. Case is
%   ignored for property names. [] is a valid OPTIONS argument.
%   
%   VAL = GNI_GET(OPTIONS,'NAME',DEFAULT) extracts the named property as
%   above, but returns VAL = DEFAULT if the named property is not specified
%   in OPTIONS. For example
%   
%       val = gni_get(options,'StepSize',1e-2);
%   
%   returns val = 1e-2 if the StepSize property is not specified in options.
%   

if nargin < 2
    error('Not enough input arguments.');
end

if ~isempty(options) && ~isa(options,'struct')
    error('First argument must be an options structure created with GNISET.');
end

if isempty(options)
  if nargin == 3
    out = default;
  else
    out = [];
  end
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
j = strncmpi(Names,name,length(name));

if max(j) == 0
      error(['Unrecognized property name ''%s''.  ' ...
                 'See GNISET for possibilities.'], name);
elseif length(find(j)) > 1            % if more than one match
  % Check for any exact matches (in case any names are subsets of others)
  k = strcmpi(name,Names);
  if max(k) == 1
    j = find(k);
  else
      msg = sprintf('Ambiguous property name ''%s'' ', name);
      l = find(j);
      msg = [msg '(' strjoin(Names(l(1),:))];
      for k = 2:length(l)
          msg = [msg ', ' strjoin(Names(l(k),:))];
      end
      msg = sprintf('%s).', msg);
      error(msg);
  end             
end

eval(['out = options.' strjoin(Names(j,:)) ';'])

if (nargin == 3) && isempty(out)
  out = default;
end