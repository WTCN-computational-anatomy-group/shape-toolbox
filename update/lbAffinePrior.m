function model = lbAffinePrior(model, opt)
% FORMAT model = lbAffinePrior(model, opt)
% model - Model-specific data
% opt   - Options
%
% All inputs are structures that can either be in memory or on disk in the
% form of a mat file. In the latter case, it is read and, if needed,
% written back.
%--------------------------------------------------------------------------
% Compute KL-divergence for the affine Wishart prior
%__________________________________________________________________________
% Copyright (C) 2018 Wellcome Centre for Human Neuroimaging

    % =====================================================================
    % Read input from disk (if needed)
    [model, modelpath, opt] = fileToStruct(model, opt);
        
    % =====================================================================
    % Compute lower bound 
    if isfinite(opt.q.n0) && opt.q.n0 > 0
        % if n0 >  0 : Bayesian update -> KL divergence
        model.lb.Aq.val = spm_prob('Wishart', 'kl', ...
                                   model.q.n, model.q.A, ...
                                   opt.q.A0,  opt.q.n0);
        model.lb.Aq.type = 'kl';
        model.lb.Aq.name = '-KL Affine precision';
    end
    % if n0 == 0   : Maximum Likelihood
    % If n0 == inf : Fixed value
                   
    % =====================================================================
    % Exit
    model = structToFile(model, modelpath);
end