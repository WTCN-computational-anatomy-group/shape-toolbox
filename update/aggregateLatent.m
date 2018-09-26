function model = aggregateLatent(dat, model, opt)
% FORMAT model = aggregateLatent(dat, model, opt)
% dat   - Subject-specific data
% model - Model-specific data
% opt   - Options
%
% All inputs are structures that can either be in memory or on disk in the
% form of a mat file. In the latter case, it is read and, if needed,
% written back.
%--------------------------------------------------------------------------
% Aggregate latent-related data:
%   >> z.z / z.zz / z.S / z.n / lb.z
%__________________________________________________________________________
% Copyright (C) 2018 Wellcome Centre for Human Neuroimaging

    % =====================================================================
    % Read input from disk (if needed)
    [dat, ~, model, modelpath, opt] = fileToStruct(dat, model, opt);

    % =====================================================================
    % Aggregate data
    model.z.z      = 0;
    model.z.zz     = 0;
    model.z.S      = 0;
    model.z.n      = 0;
    model.z.Z      = zeros(opt.pg.K, opt.v.N+opt.f.N);
    model.lb.z.val = 0;
    model.lb.z.type = 'kl';
    model.lb.z.name = '-KL Latent';
    for n=1:numel(dat)
        if iscell(dat)
            dat1 = dat{n};
        else
            dat1 = dat(n);
        end
        if isstring(dat1)
            dat1 = load(dat1);
        end
        model.z.Z(:,n) = dat1.z.z;
        model.z.n      = model.z.n      + 1;
        model.z.z      = model.z.z      + dat1.z.z;
        model.z.zz     = model.z.zz     + dat1.z.zz;
        model.z.S      = model.z.S      + dat1.z.S;
        model.lb.z.val = model.lb.z.val + dat1.z.lb.val;
    end  
        
    % =====================================================================
    % Exit
    model = structToFile(model, modelpath);
end