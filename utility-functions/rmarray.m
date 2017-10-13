function a = rmarray(a)
% FORMAT a = rmarray(a)
% Clear or remove from disk a file, file_array or Matlab array.

    if isa(a, 'file_array')
        delete(a.fname)
    elseif ischar(a)
        delete(a)
    elseif isnumeric(a)
        a = [];
    end
    
end