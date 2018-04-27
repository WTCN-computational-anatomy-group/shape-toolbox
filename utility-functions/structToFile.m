function varargout = structToFile(varargin)
% FORMAT [output1,  ...] = structToFile(struct1, path1, ...)
%
% input1  - Structur or path to a mat file
% struct1 - Structure
% path1   - Path to a mat file (or empty)
%
% Transform a mat file into a structure if needed.

    varargout = cell(1,numel(varargin)/2);
    j = 0;
    for i=1:2:numel(varargin)
        if ~isempty(varargin{i+1})
            varargout{j+1} = varargin{i+1};
            dat = varargin{i};
            save(varargin{i+1}, 'dat', '-struct');
        else
            varargout{j+1} = varargin{i};
        end
        j = j + 1;
    end

end