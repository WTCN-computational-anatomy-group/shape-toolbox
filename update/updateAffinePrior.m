function model = updateAffinePrior(model, opt)
% FORMAT dat = updateAffinePrior(dat, model, opt)
% dat   - Subject-specific data
% model - Model-specific data
% opt   - Options
%
% All inputs are structures that can either be in memory or on disk in the
% form of a mat file. In the latter case, it is read and, if needed,
% written back.
%--------------------------------------------------------------------------
% Update affine Wishart prior
%__________________________________________________________________________
% Copyright (C) 2018 Wellcome Centre for Human Neuroimaging

    % =====================================================================
    % Read input from disk (if needed)
    [model, modelpath, opt] = fileToStruct(model, opt);
        
    % =====================================================================
    % Update 
    if opt.q.n0 == inf
        % If n0 == inf : Fixed value
        model.q.A       = opt.q.A0;
        model.q.n       = inf;
        model.q.LogDetA = spm_matcomp('LogDet', model.q.A);
    elseif opt.q.n0 == 0
        % if n0 == 0 : Maximum Likelihood
        rind            = opt.q.rind;
        model.q.A       = spm_matcomp('Inv', (model.q.qq(rind,rind) + model.q.S(rind,rind))/model.q.n);
        model.q.LogDetA = spm_matcomp('LogDet', model.q.A);
    else
        % if n0 >  0 : Bayesian update
        rind            = opt.q.rind;
        model.q.A       = spm_prob('Wishart', 'up', ...
            model.q.n, 0, model.q.qq(rind,rind) + model.q.S(rind,rind), ...
            opt.q.A0, opt.q.n0);
        model.q.n       = model.q.n + opt.q.n0;
        model.q.LogDetA = spm_prob('Wishart', 'ELogDet', model.q.A, model.q.n);
    end
                   
    % =====================================================================
    % Exit
    model = structToFile(model, modelpath);
end