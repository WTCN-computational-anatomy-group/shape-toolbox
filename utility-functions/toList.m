function list = toList(s, fields, varargin)
%__________________________________________________________________________
% Extract a value from all elements of a struct array or a cell array of 
% structures and return it in a cell array.
%--------------------------------------------------------------------------
% FORMAT list = toList(s, fields, ...)
%
% REQUIRED
% --------
% s       - Structure array or cell array of structures
% fields  - String describing the value to extract
%
% KEYWORD ARGUMENTS
% -----------------
% missing  - value to use if the target field does not exist [zeros(0)]
% function - Function to apply to each extracted value [identity]
%--------------------------------------------------------------------------
%   Ex:
%       list = toList(s, '.dat.vs(1)', 'missing', 1)
%__________________________________________________________________________

    % --- Parse inputs
    p = inputParser;
    p.FunctionName = 'toList';
    p.addRequired('s',      @(X) iscell(X) || isstruct(X));
    p.addRequired('fields', @ischar);
    p.addParameter('missing', zeros(0));
    p.addParameter('function',    @(X) X, @(X) isa(X, 'function_handle'));
    p.parse(s, fields, varargin{:});

    N = numel(s);
    list = cell(size(s));
    if iscell(s)
        for n=1:N
            list{n} = p.Results.function(defval(s{n}, fields, p.Results.missing));
        end
    elseif isstruct(s)
        for n=1:N
            list{n} = p.Results.function(defval(s(n), fields, p.Results.missing));
        end
    end

end