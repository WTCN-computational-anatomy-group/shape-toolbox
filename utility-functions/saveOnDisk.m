function out = saveOnDisk(where, obj, varargin)
% FORMAT out = createDiskOutput(where, obj,
%                               ('append'),
%                               ('name', name),
%                               ('mat', mat),
%                               ('mat0', mat0),
%                               ('type', type))
% where    - file_array or filename (.nii or .mat)
% obj      - Object to save
% 'append' - If output is a .mat and 'append' is set, appends the object
%            to the file [false]
% name     - If output is a .mat, variable name ['obj']
% mat      - If output is .nii, Affine matrix stored in the header [eye(4)]
% type     - If output is .nii, on-disk data type ['float32']
% out      - Output object (if .nii -> the file_array, else -> the object)
%
% Save an object on disk, whether it as a Nifti image file or as a Matlab 
% mat file.

    p = inputParser;
    p.addRequired('where');
    p.addRequired('obj');
    p.addParameter('name', 'obj'); % Unused now, need to remove it
    p.addParameter('mat', eye(4));
    p.addParameter('mat0', eye(4));
    p.addParameter('type', 'float32');
    p.parse(where, obj, varargin{:});
    
    if ischar(where)
        if endsWith(where, '.nii')
            out = file_array(where);
            out.dtype = p.Results.type;
            n = nifti;
            n.dat = out;
            n.mat = p.Results.mat;
            n.mat0 = p.Results.mat0;
            create(n);
            out = n.dat;
        else
            error('Extension not handled. Must be .nii or .mat')
        end
    elseif isa(where, 'file_array')
        out = where;
    else
        out = obj;
        return
    end
    
    out.dim = size(obj);
    out(:) = obj(:);
    
end