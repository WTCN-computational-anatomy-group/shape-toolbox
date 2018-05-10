function model = aggregateVelocity(dat, model, opt)
% FORMAT model = aggregateVelocity(dat, model, opt)
% dat   - Subject-specific data
% model - Model-specific data
% opt   - Options
%
% All inputs are structures that can either be in memory or on disk in the
% form of a mat file. In the latter case, it is read and, if needed,
% written back.
%--------------------------------------------------------------------------
% Aggregate latent-related data:
%   >> v.n / v.tr / v.reg / v.uncty / lb.v1 / lb.v2 
%__________________________________________________________________________
% Copyright (C) 2018 Wellcome Centre for Human Neuroimaging

    % =====================================================================
    % Read input from disk (if needed)
    [dat, ~, model, modelpath, opt] = fileToStruct(dat, model, opt);

    % =====================================================================
    % Aggregate data
    nf = 0;
    nv = 0;
    model.v.n       = 0; % For lambda update
    model.v.tr      = 0; % For lambda update
    model.v.reg     = 0; % For lambda update
    model.v.uncty   = 0; % For lambda update
    if opt.v.N
        model.lb.v1.val = 0; % Log-likelihood (observed velocity)
        model.lb.v1.type = 'll';
        model.lb.v1.name = 'LL Velocity';
    end
    if opt.f.N
        model.lb.v2.val = 0; % KL-divergence  (latent velocity)
        model.lb.v2.type = 'kl';
        model.lb.v2.name = '-KL Velocity';
    end
    for n=1:numel(dat)
        if iscell(dat)
            dat1 = dat{n};
        else
            dat1 = dat(n);
        end
        if isstring(dat1)
            dat1 = load(dat1);
        end
        if defval(dat1.f, '.observed', true)
            nf = nf + 1;
            model.v.tr    = model.v.tr    + dat1.v.lb.tr;
            model.v.uncty = model.v.uncty + dat1.v.lb.uncty;
            model.lb.v2.val = model.lb.v2.val + dat1.v.lb.val;
        else
            nv = nv + 1;
            model.lb.v1.val = model.lb.v1.val + dat1.v.lb.val;
        end
        model.v.n      = model.v.n   + 1;
        model.v.reg    = model.v.reg + dat1.v.lb.reg;
    end
        
    % =====================================================================
    % Exit
    model = structToFile(model, modelpath);
end