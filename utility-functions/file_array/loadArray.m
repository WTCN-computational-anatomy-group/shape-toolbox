function a = loadArray(h, v1, v2)
% FORMAT a = loadArray(fa, ['copy', dtype])
% h      - Nifti / HandleArray / file_array / array
% 'copy' - If specified and the input is in memory, forces a full copy.
% dtype  - Force the output array to be of class |dtype|.
% a      - In-memory HandleArray. Data is in a.dat
%
% Load in a HandleArray the data contained in a file_array.
% If the input is already in memory, just return the input.
    
    % --- Parse arguments
    switch nargin
        case 3
            if strcmp(v1, 'copy')
                force_copy = true;
                dtype = v2;
            else
                dtype = v1;
                force_copy = strcmp(v2, 'copy');
            end
        case 2
            if strcmp(v1, 'copy')
                force_copy = true;
                dtype = 0;
            else
                dtype = v1;
                force_copy = false;
            end
        case 1
            force_copy = false;
            dtype = 0;
    end
    
    if isa(h.dat, 'file_array')
    %--- Data on disk
        a = HandleArray;
        if issame(dtype, 0)
            dtype = fileArray2MatlabType(h.dat.dtype);
        end
        a.dat = zeros(size(h.dat), dtype);
        a.dat(:) = h.dat(:);
    else
    %--- Data in memory
        if issame(dtype, 0)
            dtype = class(a.dat);
        end
        if force_copy || dtype ~= class(a.dat)
            a = HandleArray;
            a.dat = cast(fa.dat, dtype);
        else
            a = h;
        end
    end
end