function model = aggregateAffine(dat, model, ~)
% FORMAT model = aggregateAffine(dat, model, opt)
% dat   - Subject-specific data
% model - Model-specific data
% opt   - Options
%
% All inputs are structures that can either be in memory or on disk in the
% form of a mat file. In the latter case, it is read and, if needed,
% written back.
%--------------------------------------------------------------------------
% Aggregate affine-related data:
%   >> q.q / q.qq / q.S / q.n / lb.q
%__________________________________________________________________________
% Copyright (C) 2018 Wellcome Centre for Human Neuroimaging

    % =====================================================================
    % Read input from disk (if needed)
    [dat, ~, model, modelpath] = fileToStruct(dat, model);

    % =====================================================================
    % Aggregate data
    model.q.q      = 0;
    model.q.qq     = 0;
    model.q.S      = 0;
    model.q.n      = 0;
    model.lb.q.val = 0;
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
            model.q.n      = model.q.n      + 1;
            model.q.q      = model.q.q      + dat1.q.q;
            model.q.qq     = model.q.qq     + dat1.q.qq;
            model.q.S      = model.q.S      + dat1.q.S;
            model.lb.q.val = model.lb.q.val + model.q.lb.val;
        end
    end  
        
    % =====================================================================
    % Exit
    model = structToFile(model, modelpath);
end