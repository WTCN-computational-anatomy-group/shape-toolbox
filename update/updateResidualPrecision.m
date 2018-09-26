function model = updateResidualPrecision(model, opt)
% FORMAT dat = updateResidualPrecision(dat, model, opt)
% dat   - Subject-specific data
% model - Model-specific data
% opt   - Options
%
% All inputs are structures that can either be in memory or on disk in the
% form of a mat file. In the latter case, it is read and, if needed,
% written back.
%--------------------------------------------------------------------------
% Update residual precision
%__________________________________________________________________________
% Copyright (C) 2018 Wellcome Centre for Human Neuroimaging

    % =====================================================================
    % Read input from disk (if needed)
    [model, modelpath, opt] = fileToStruct(model, opt);
        
    % =====================================================================
    % Update 
    if opt.v.n0 == inf
        % If n0 == inf : Fixed value
        model.v.l       = opt.v.l0;
        model.v.n       = inf;
        model.v.LogLambda = log(model.v.l);
    elseif opt.v.n0 == 0
        % if n0 == 0 : Maximum Likelihood
        K = 3*prod(opt.tpl.lat);
        model.v.n = model.mixreg.w(1) * (opt.v.N+opt.f.N);
        model.v.l = (model.mixreg.w(1)/K)*(model.v.uncty + model.v.tr + model.v.reg);
        model.v.l = model.v.n/model.v.l;
        model.v.LogLambda = log(model.v.l);
    else
        % if n0 >  0 : Bayesian update
        K = 3*prod(opt.tpl.lat);
        model.v.n = opt.v.n0 + model.mixreg.w(1) * (opt.v.N+opt.f.N);
        model.v.l = opt.v.n0/opt.v.l0 + (model.mixreg.w(1)/K)*(model.v.uncty + model.v.tr + model.v.reg);
        model.v.l = model.v.n/model.v.l;
        model.v.LogLambda = spm_prob('Gamma', 'ELog', model.v.l, model.v.n, K, 'normal');
    end
                   
    % =====================================================================
    % Exit
    model = structToFile(model, modelpath);
end