%... The MatMol Group (2016)
    function varargout = funcclapeyron(varargin)
%... Implementation of Clapeyron equation
    
    K = 8314/(2791.2e3*18);

%... Depending on the output variable, different versions of
%... the formula
    switch varargin{2}
        case {'T'}
            P = varargin{1};
            T = 1./( 1/273.11 - K*log(P/611.73));
            varargout{1} = T;
        case {'P'}
            T = varargin{1};
            P = 611.73*exp(1/K*(1/273.11-1./T));
            varargout{1} = P;
    end