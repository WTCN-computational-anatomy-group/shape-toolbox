function deallocateArray(a, rm)
% FORMAT deallocateArray(a)
% a     - Nifti or HandleArray object.
% (rm)  - Remove file if storage is on disk (default:false)
% 
% * If data is in memory, replace the array with an empty one (keep dtype
%   the same).
% * If |rm| is true and data is on disk, remove the file.
% * In all cases, handle object are not deleted (nifti, HandleArray, 
%   file_array)
    if nargin < 2
        rm = false;
    end
    if ~isa(a.dat, 'file_array')
        a.dat = zeros(0, class(a.dat));
    elseif rm
        delete a.dat.fname
    end
end