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
    end

    % =====================================================================
    % Latent coordinates
    if opt.optimise.z.z
        dat = initLatent(dat, model, opt);
        dat = lbLatent(dat, model, opt);
    end

    % =====================================================================
    % Velocity / transforms
    if opt.optimise.v.v
        dat = initVelocity(dat, model, opt);
    end
    
    % ====================================================================
    % Tissue scales
    dat = initTissueScale(dat, model, opt);
    
    % ====================================================================
    % Template
    if opt.optimise.tpl.a
        switch lower(opt.tpl.update)
            case 'ml'
                dat = initPush(dat, model, opt);
            case 'map'
                dat = updateWarpedTemplate(dat, model, opt);
                dat = updateTemplateGradHess(dat, model, opt);
        end
    end
    
    % =====================================================================
    % Exit
    dat = structToFile(dat, datpath);
end
