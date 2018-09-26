function dat = initVelocity(dat, ~, opt)
% FORMAT dat = initVelocity(dat, model, opt)
% dat   - Subject-specific data
% model - Model-specific data
% opt   - Options
%
% All inputs are structures that can either be in memory or on disk in the
% form of a mat file. In the latter case, it is read and, if needed,
% written back.
%--------------------------------------------------------------------------
% Update affine parameters by Gauss-Newton
%__________________________________________________________________________
% Copyright (C) 2018 Wellcome Centre for Human Neuroimaging

    % =====================================================================
    % Read input from disk (if needed)
    [dat, datpath, opt] = fileToStruct(dat, opt);

    % =====================================================================
    % If the velocity is an observed -> nothing to do
    if defval(dat.v, '.observed', false)
        dat = structToFile(dat, datpath);
        return
    end
    
    % =====================================================================
    % Velocity (zero)
    dat.v.v   = prepareOnDisk(dat.v.v, [opt.tpl.lat 3]);
    dat.v.v(:) = 0;
    dat.v.ok  = 1;
    dat.v.ok2 = 1;
    
    
    % =====================================================================
    % Diffeomorphic transform (identity)
    id = spm_warps('identity', opt.tpl.lat);
    dat.v.iphi = saveOnDisk(dat.v.iphi, id);
    
    % =====================================================================
    % Complete transform (compose)
    if isfield(dat, 'q') && isfield(dat.q, 'A'),  A = dat.q.A;
    else,                                         A = eye(4);  end
    dimf   = [size(dat.f.f) 1 1];
    latf   = dimf(1:3);
    dat.v.ipsi = reconstructIPsi(...
        A, ...                      % Affine transform
        id, ...                     % Diffeomorphic transform
        'lat',    latf, ...         % Image lattice
        'Mf',     dat.f.M, ...      % Image voxel-to-world
        'Mmu',    opt.tpl.M, ...    % Template lattice
        'output', dat.v.ipsi, ...   % Output file_array (save disk)
        'debug',  opt.ui.debug);    % Print debugging stuff. (usually no)
    clear id
              
    % =====================================================================
    % Exit
    dat = structToFile(dat, datpath);
end