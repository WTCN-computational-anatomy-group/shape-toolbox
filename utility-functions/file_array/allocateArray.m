function allocateArray(a, dim)
% FORMAT a = allocateArray(a, dim)
% a     - Nifti or HandleArray (data is in a.dat).
% (dim) - New dimension of the array. (default: keep existing)
%
% (Re)allocate array, which can be on disk or in memory.
% This function handles a referenced array: the input is *always* modified.
    if nargin < 2
        dim = size(a.dat);
    end
    if isa(a.dat, 'file_array')
    % Data on disk
        if ~issame(size(a.dat), dim)
            a.dat.dim = dim;
        end
        % If the array is not allocated (enough), allocate
        try
            a.dat(end);
        catch
            a.dat(:) = zeros(numel(a.dat), a.dat.dtype);
        end
    else
    % Data in memory
        if ~issame(numel(a.dat), prod(dim))
        % Total number of element is different: reallocate
            dtype = class(a.dat);
            a.dat = zeros(dim, dtype);
        elseif ~issame(size(a.dat), dim)
        % Same number of element but different shape: reshape
            a.dat = reshape(a.dat, dim);
        end
    end
end