function dat = lbLatent(dat, model, opt)
% FORMAT dat = lbLatent(dat, model, opt)
% dat   - Subject-specific data
% model - Model-specific data
% opt   - Options
%
% All inputs are structures that can either be in memory or on disk in the
% form of a mat file. In the latter case, it is read and, if needed,
% written back.
%--------------------------------------------------------------------------
% Update latent coordinates KL divergence (dat.z.lb.val)
%__________________________________________________________________________
% Copyright (C) 2018 Wellcome Centre for Human Neuroimaging

    % =================================================================
    % Read input from disk (if needed)
    [dat, datpath, model, ~, opt] = fileToStruct(dat, model, opt);

    % =================================================================
    % KL-divergence
    dat.z.lb.val  = -0.5*( trace((dat.z.S + dat.z.zz) * model.z.A) ...
                           - model.z.LogDetA ...
                           - spm_matcomp('LogDet', dat.z.S) ...
                           - opt.pg.K );
    dat.z.lb.type = 'kl';
    dat.z.lb.name = '-KL Latent';

    % =================================================================
    % Exit
    dat = structToFile(dat, datpath);
        
end