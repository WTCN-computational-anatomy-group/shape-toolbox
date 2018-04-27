function dat = shape_init_subject(dat, model, opt)
% FORMAT dat = shape_init_subject(dat, model, opt)
% dat   - Subject-specific data
% model - Model-specific data
% opt   - Options
%
% All inputs are structures that can either be in memory or on disk in the
% form of a mat file. In the latter case, it is read and, if needed,
% written back.
%--------------------------------------------------------------------------
% Init subject-specific variables.
% These parameters cannot be provided and must be optimised.
% This function can be run as an independent job on a cluster.
%__________________________________________________________________________
% Copyright (C) 2018 Wellcome Centre for Human Neuroimaging

    % =====================================================================
    % Read input from disk (if needed)
    [dat, datpath, model, ~, opt] = fileToStruct(dat, model, opt);
    
    % =====================================================================
    % Affine coordinates
    if defval(opt.f, '.observed', true) && opt.optimise.q.q
        dat = initAffine(dat, model, opt);
        dat = lbAffine(dat, model, opt);
    end

    % =====================================================================
    % Latent coordinates
    dat = initLatent(dat, model, opt);
    dat = lbLatent(dat, model, opt);

    % =====================================================================
    % Velocity / transforms
    dat = initVelocity(dat, model, opt);
    dat = lbVelocity(dat, model, opt);
    
    % ====================================================================
    % Template
    switch lower(opt.tpl.update)
        case 'ml'
            dat = initPush(dat, model, opt);
        case 'map'
            dat = initTemplateGradHess(dat, model, opt);
            dat = lbMatching(dat, model, opt);
    end
    
    % =====================================================================
    % Exit
    dat = structToFile(dat, datpath);
end
