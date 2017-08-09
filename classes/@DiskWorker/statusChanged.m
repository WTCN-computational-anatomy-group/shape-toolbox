function statusChanged(obj, varargin)
% FORMAT obj.statusChanged('prop1', 'prop2', ...)
% 
% Sets all downstream dependencies of the input properties to utd == false 
% and all the input properties to utd == true.

    if nargin <= 1
        return
    end
    for i=1:numel(varargin)
        if isfield(obj.dependencies, varargin{i})
            if isfield(obj.downdep, varargin{i})
                downstream = obj.downdep.(varargin{i});
            else
                downstream = {};
            end
            obj.statusChanged(downstream{:});
            for j=1:numel(downstream)
                if isfield(obj.utd, downstream{j})
                    obj.utd.(downstream{j}) = false;
                end
            end
        end
    end
    for i=1:numel(varargin)
        obj.utd.(varargin{i}) = true;
    end
end