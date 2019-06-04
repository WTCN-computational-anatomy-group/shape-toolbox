function model = lbResidualPrior(model, opt)
% FORMAT model = lbResidualPrior(model, opt)
% model - Model-specific data
% opt   - Options
%
% All inputs are structures that can either be in memory or on disk in the
% form of a mat file. In the latter case, it is read and, if needed,
% written back.
%--------------------------------------------------------------------------
% Compute KL-divergence for the latent residual precision
%__________________________________________________________________________
% Copyright (C) 2018 Wellcome Centre for Human Neuroimaging

    % =====================================================================
    % Read input from disk (if needed)
    [model, modelpath, opt] = fileToStruct(model, opt);
        
    % =====================================================================
    % Compute lower bound 
    if isfinite(opt.v.n0) && opt.v.n0 > 0
        % if n0 >  0 : Bayesian update -> KL divergence
        model.lb.l.val  = -spm_prob('Gamma', 'kl', ...
                                    model.v.l, model.v.n, ...
                                    opt.v.l0,  opt.v.n0, ...
                                    prod(opt.tpl.lat)*3, 'normal');
        model.lb.l.type = 'kl';
        model.lb.l.name = '-KL Residual precision';
    end
    % if n0 == 0   : Maximum Likelihood
    % If n0 == inf : Fixed value
                   
    % =====================================================================
    % Exit
    model = structToFile(model, modelpath);
end