function array = toArray(s, fields, varargin)
%__________________________________________________________________________
% Extract a value from all elements of a struct array or a cell array of 
% structures and return it in an array.
%--------------------------------------------------------------------------
% FORMAT list = toArray(s, fields, ...)
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
%       array = toArray(s, '.dat.vs(1)', 'missing', 1)
%__________________________________________________________________________

    list = toList(s, fields, varargin{:});
    array = cell2mat(list);

end