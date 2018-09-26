function dat = shape_init_subject2(dat, model, opt)
% FORMAT dat = shape_init_subject2(dat, model, opt)
% dat   - Subject-specific data
% model - Model-specific data
% opt   - Options
%
% All inputs are structures that can either be in memory or on disk in the
% form of a mat file. In the latter case, it is read and, if needed,
% written back.
%--------------------------------------------------------------------------
% Init subject-specific variables that **require a template**.
% These parameters cannot be provided and must be optimised.
% This function can be run as an independent job on a cluster.
%__________________________________________________________________________
% Copyright (C) 2018 Wellcome Centre for Human Neuroimaging

    % =====================================================================
    % Read input from disk (if needed)
    [dat, datpath, model, ~, opt] = fileToStruct(dat, model, opt);
    
    
    % =====================================================================
    % Pull template
    if defval(opt.f, '.observed', true)
        dat = updateWarpedTemplate(dat, model, opt); % Pull and reconstruct
        dat = lbMatching(dat, model, opt);           % Update matching term
    end
    
    % =====================================================================
    % Affine coordinates
    if defval(opt.f, '.observed', true) && opt.optimise.q.q
        dat = initAffineLaplace(dat, model, opt); % Hessian for Laplace approximation
        dat = lbAffine(dat, model, opt);
    end

    % =====================================================================
    % Velocity / transforms
    if defval(opt.f, '.observed', true)
        dat = initVelocityShapeLaplace(dat, model, opt); % Hessian for Laplace approximation
        dat = lbVelocityShape(dat, model, opt);
    end
    
    % =====================================================================
    % Exit
    dat = structToFile(dat, datpath);
end
