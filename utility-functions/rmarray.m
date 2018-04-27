function a = rmarray(a)
% FORMAT a = rmarray(a)
% Clear or remove from disk a file, file_array or Matlab array.

    if isa(a, 'file_array')
        names = {a.fname};
        for i=1:numel(names)
            if exist(names{i}, 'file')
                delete(names{i})
            end
        end
    elseif ischar(a)
        if exist(a, 'file')
            delete(a)
        end
    elseif iscell(a)
        for i=1:numel(a)
            a{i} = rmarray(a{i});
        end
    elseif isnumeric(a)
        a = [];
    end
    
end