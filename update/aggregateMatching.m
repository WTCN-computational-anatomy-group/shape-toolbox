function model = aggregateMatching(dat, model, ~)
% FORMAT model = aggregateMatching(dat, model, opt)
% dat   - Subject-specific data
% model - Model-specific data
% opt   - Options
%
% All inputs are structures that can either be in memory or on disk in the
% form of a mat file. In the latter case, it is read and, if needed,
% written back.
%--------------------------------------------------------------------------
% Aggregate matching terms:
%   >> lb.f
%__________________________________________________________________________
% Copyright (C) 2018 Wellcome Centre for Human Neuroimaging

    % =====================================================================
    % Read input from disk (if needed)
    [dat, ~, model, modelpath] = fileToStruct(dat, model);

    % =====================================================================
    % Aggregate data
    model.lb.f.val = 0;
    model.lb.f.type = 'll';
    model.lb.f.name = 'Matching likelihood';
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
            model.lb.f.val = model.lb.f.val + dat1.f.lb.val;
        end
    end  
        
    % =====================================================================
    % Exit
    model = structToFile(model, modelpath);
end