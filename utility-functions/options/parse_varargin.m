function out = parse_varargin(in, out)
% FORMAT opt = parse_varargin(in, opt)
% in    - varargin cell
% (opt) - (Optional) Structure of options
%
% This helper can be called at the beginning of a function that takes
% variable arguments.
% - The varargin variable should be passed as input and must be of format: 
%   {'key1', value1, 'key2', value2, ...}. It is parsed and used to build a 
%   structure associating keys with values.
% - If a key is of the form 'key.otherkey', the hierarchy is kept in the 
%   output structure s.t. opt.key.otherkey = value.
% - If an initial opt structure is given as an input, it is only updated 
%   with varargin values (i.e., keys that are not specified in the varargin 
%   object are not changed). 
% - If the number of cells in the varargin object is odd, and if the last
%   one is a struct, fields from this struct are also used to update the
%   input structure.
%
% USAGE
% varargin = {'position', [3 3]};
% dft = struct('color', 'red', 'position', [1 1]); % Default values
% opt = parse_varargin(varargin, dft); % Parse variable arguments
% dft
% opt
%
% > dft = 
% > 
% >   struct with fields:
% > 
% >        color: 'red'
% >     position: [1 1]
% > opt = 
% > 
% >   struct with fields:
% > 
% >        color: 'red'
% >     position: [3 3]
    if nargin < 2
        out = struct;
    end
    for i=1:floor(numel(in)/2)
        key = in{i*2-1};
        value = in{i*2};
        % In case of option tree (i.e., struct of struct of ...)
        key = strsplit(key, '.');
        for j=1:length(key)
            key{j} = char(key{j});
        end
        for j=1:length(key)-1
            if j == 1
                prev = out;
            else
                prev = out.(key{1:j-1});
            end
            if ~isfield(prev, key{j}) || ~isstruct(prev.(key{j}))
                out.(key{1:j}) = struct;
            end
        end
        out = setfield(out, key{:}, value);
    end
    % If odd number of arguments, the last one can be an option structure
    if mod(numel(in), 2) == 1 && isstruct(in{end})
        out = update_struct(in{end}, out);
    end
end