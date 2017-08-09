function n = checkNifti(n, dim)
% FORMAT checkNifti(n, (dim))
% n     - Nifti object
% dim   - Dimensions to set (optional)
%
% Check that the nifti file is ready to be written in:
% 1) If the file does not exists, call create(n)
% 2) If needed, update the array dimensions
% 3) Just incase, initialise(n.dat)

    if nargin > 1
        n.dat.dim = dim;
    end
    if ~exist(n.dat.fname, 'file')
        create(n);
    end
    initialise(n.dat);
    
end