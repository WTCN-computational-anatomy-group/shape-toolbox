function model = updateLatentPrior(model, opt)
% FORMAT dat = updateLatentPrior(dat, model, opt)
% dat   - Subject-specific data
% model - Model-specific data
% opt   - Options
%
% All inputs are structures that can either be in memory or on disk in the
% form of a mat file. In the latter case, it is read and, if needed,
% written back.
%--------------------------------------------------------------------------
% Update latent Wishart prior
%__________________________________________________________________________
% Copyright (C) 2018 Wellcome Centre for Human Neuroimaging

    % =====================================================================
    % Read input from disk (if needed)
    [model, modelpath, opt] = fileToStruct(model, opt);
        
    % =====================================================================
    % Update 
    if opt.z.n0 == inf
        % If n0 == inf : Fixed value
        model.z.A       = opt.z.A0;
        model.z.n       = inf;
        model.z.LogDetA = spm_matcomp('LogDet', model.z.A);
    elseif opt.z.n0 == 0
        % if n0 == 0 : Maximum Likelihood
        model.z.A       = spm_matcomp('Inv', (model.z.zz + model.z.S)/model.z.n);
        model.z.LogDetA = spm_matcomp('LogDet', model.z.A);
    else
        % if n0 >  0 : Bayesian update
        model.z.A       = spm_prob('Wishart', 'up', ...
            model.z.n, 0, model.z.zz + model.z.S, ...
            opt.z.A0, opt.z.n0);
        model.z.n       = model.z.n + opt.z.n0;
        model.z.LogDetA = spm_prob('Wishart', 'ELogDet', model.z.A, model.z.n);
    end
                   
    % =====================================================================
    % Exit
    model = structToFile(model, modelpath);
end