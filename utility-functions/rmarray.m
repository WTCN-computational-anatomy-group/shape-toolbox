function a = rmarray(a)
% FORMAT a = rmarray(a)
% Clear or remove from disk a file, file_array or Matlab array.

    if isa(a, 'file_array')
        if exist(a.fname, 'file')
            delete(a.fname)
        end
    elseif ischar(a)
        if exist(a, 'file')
            delete(a)
        end
    elseif isnumeric(a)
        a = [];
    end
    
end