function model = lbSubspace(model, opt)
% FORMAT model = lbSubspace(model, opt)
% model - Model-specific data
% opt   - Options
%
% All inputs are structures that can either be in memory or on disk in the
% form of a mat file. In the latter case, it is read and, if needed,
% written back.
%--------------------------------------------------------------------------
% Update prior log-likelihood for the principal subspace
%__________________________________________________________________________
% Copyright (C) 2018 Wellcome Centre for Human Neuroimaging

    % =====================================================================
    % Read input from disk (if needed)
    [model, modelpath, opt] = fileToStruct(model, opt);
     
    % =====================================================================
    % Lower Bound
    model.lb.w.val = llPriorSubpsace(model.pg.w, model.pg.ww, opt.pg.LogDetL);
    model.lb.w.type = 'll';
    model.lb.w.name = 'Subspace prior';
                   
    % =====================================================================
    % Exit
    model = structToFile(model, modelpath);
end