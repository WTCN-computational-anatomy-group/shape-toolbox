function varargout = fileToStruct(varargin)
% FORMAT [struct1, path1, ...] = fileToStruct(input1, ...)
%
% input1  - Structur or path to a mat file
% struct1 - Structure
% path1   - Path to a mat file (or empty)
%
% Transform a mat file into a structure if needed.

    varargout = cell(1,2*numel(varargin));
    j = 0;
    for i=1:numel(varargin)
        if ischar(varargin{i})
            varargout{j+2} = varargin{i};
            varargout{j+1} = load(varargin{i});
        else
            varargout{j+2} = '';
            varargout{j+1} = varargin{i};
        end
        j = j + 2;
    end

end