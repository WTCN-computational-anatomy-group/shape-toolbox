function a = initArray(dim, dtype, fname, opt)
% FORMAT a = initArray(dim, dtype, storage, nifti, opt)
% dim     - Dimensions to allocate. Can be 0 if the allocation is
%           delayed to a later call.
% dtype   - Array data type ('single', 'double', 'uint8t', ...).
% fname   - 'memory.nii': data is stored in RAM as a (virtual) nifti
%           'memory*':    data is stored in RAM as a HandleArray(array)
%           'fname.nii':  data is stored on disk as a nifti
%           'fname.*':    data is stored on disk as a
%                         HandleArray(file_array)
% opt     - Structure of nifti options.
%
% Initialize the array depending on its memory type (file or
% memory, nifti or simple file_array).
%
% If the 
    if nargin < 5
        opt = struct;
    end
    if endsWith(fname, '.nii')
        a = nifti;
        if ~startsWith(fname, 'memory')
            a.dat       = file_array;
            a.dat.fname = fname;
            a.dat.dtype = dtype;
            if ~issame(dim, 0)
                a.dat.dim = dim;
            end
        else
            a.dat = zeros(dim, dtype);
        end
        fields = fieldnames(opt);
        for i=1:length(fields)
            switch fields(i)
                case 'intent'
                    a.intent = opt.intent;
                case 'mat'
                    a.mat = opt.mat;
                case 'mat_intent'
                    a.mat_intent = opt.mat_intent;
                case 'mat0'
                    a.mat0 = opt.mat0;
                case 'mat0_intent'
                    a.mat0_intent = opt.mat0_intent;
                case 'code'
                    a.code = opt.code;
                case 'param'
                    a.param = opt.param;
            end
        end
        if ~startsWith(fname, 'memory')
            create(a);
        end
    else
        a = HandleArray;
        if ~startsWith(fname, 'memory')
            a.dat       = file_array;
            a.dat.fname = fname;
            a.dat.dtype = dtype;
            if ~issame(dim, 0)
                a.dat.dim = dim;
            end
        else
            a.dat = zeros(dim, dtype);
        end
    end
end