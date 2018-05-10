function dat = initLatent(dat, model, opt)
% FORMAT dat = initLatent(dat, model, opt)
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
    [dat, datpath, model, ~, opt] = fileToStruct(dat, model, opt);

    % =====================================================================
    % Initial values
    switch lower(opt.z.init)
        case 'rand'
            dat.z.z   = mvnrnd(zeros(opt.pg.K, 1), spm_matcomp('Inv',model.z.A));
        case 'zero'
            dat.z.z   = zero(opt.pg.K, 1);
    end
    dat.z.zz  = dat.z.z * dat.z.z';
    dat.z.S   = inv(model.v.l * model.pg.ww + model.z.A);
    
    
    % =====================================================================
    % Exit
    dat = structToFile(dat, datpath);
end