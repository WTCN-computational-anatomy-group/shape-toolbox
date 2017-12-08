function out = prepareOnDisk(where, dim, varargin)
% FORMAT out = prepareOnDisk(where, dim,
%                            ('mat', mat),
%                            ('mat0', mat0),
%                            ('type', type))
% where    - file_array or filename (.nii or .mat) or empty
% dim      - Dilensions of the object to save
% mat      - If output is .nii, Affine matrix stored in the header [eye(4)]
% type     - If output is .nii, on-disk data type ['float32']
% out      - Output object (if .nii -> the file_array, else -> the object)
%
% If where is a file_array
%   > Set the appropriate dim and return the file_array
% If where is a path.nii
%   > Prepare the appropriate Nifti file on disk and return the file_array
% If where is a path.mat or empty
%   > Allocate an appropriate Matlab array.


    p = inputParser;
    p.addRequired('where');
    p.addRequired('dim');
    p.addParameter('mat', eye(4));
    p.addParameter('mat0', eye(4));
    p.addParameter('type', 'float32');
    p.parse(where, dim, varargin{:});
    
    
    if isempty(where) || ischar(where) || isnumeric(where)
        
        if isempty(where) || isnumeric(where)
            type = p.Results.type;
            if startsWith(lower(type), 'float32')
                type = 'single';
            elseif startsWith(lower(type), 'float64')
                type = 'double';
            end
            out = zeros(dim, type);
            return
            
        elseif endsWith(where, '.nii')
            type = p.Results.type;
            if startsWith(lower(type), 'single')
                type = 'float32';
            elseif startsWith(lower(type), 'double')
                type = 'float64';
            end
            out = file_array(where);
            out.dtype = type;
            n = nifti;
            n.dat = out;
            n.mat = p.Results.mat;
            n.mat0 = p.Results.mat0;
            create(n);
            out = n.dat;
            
        else
            error('Extension not handled. Must be .nii')
            
        end
        
    elseif isa(where, 'file_array')
        out = where;
        
    else
        error('Target not handled. Must be a matlab array, file_array or path')
    end
    
    out.dim = dim;
    
end