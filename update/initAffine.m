function dat = initAffine(dat, ~, opt)
% FORMAT dat = initAffine(dat, model, opt)
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
    % Initial values (zero)
    dat.q.ok  = 1;
    dat.q.ok2 = 1;
    dat.q.q   = zeros(opt.q.M, 1);
    dat.q.qq  = zeros(opt.q.M);
    dat.q.S   = zeros(opt.q.M);
    dat.q.A   = eye(4);
              
    % =====================================================================
    % Exit
    dat = structToFile(dat, datpath);
end