function out = copyarray(in, out)
% FORMAT out = copyarray(in, (out))
% in  - array or file_array to copy
% out - array/file_array/filename where to copy
%
% Intelligent copy of (file) arrays:
% - If `in` and `out` are numeric
%       normal copy (out = in)
% - If `in` is a file_array and `out` is numeric
%       load in memory (out = numeric(in))
% - If `in` is numeric and `out` is a file_array
%       prepare file_array and write on disk (out(:) = in(:))
% - If `in` and `out` are the same file_array
%       swap them (to update dimension field, etc.)
% - If `in` and `out` are different file_arrays
%       prepare output and copy (in(:,:,z,:) = out(:,:,z,:))
%       [split copy along third dimension to save memory]

    if nargin < 2 || ~isa(out, 'file_array')
        if ~isa(in, 'file_array')
            out = in;
        else
            out = numeric(in);
        end
    else
        if isa(in, 'file_array')
            out = prepareOnDisk(out, size(in), 'type', in.dtype);
            for z=1:size(in, 3)
                out(:,:,z,:) = in(:,:,z,:);
            end
        else
            out = prepareOnDisk(out, size(in));
            out(:) = in(:);
        end
    end
    
end
