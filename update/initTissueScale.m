function dat = initTissueScale(dat, ~, opt)
% FORMAT dat = initTissueScale(dat, model, opt)
% dat   - Subject-specific data
% model - Model-specific data
% opt   - Options
%
% All inputs are structures that can either be in memory or on disk in the
% form of a mat file. In the latter case, it is read and, if needed,
% written back.
%--------------------------------------------------------------------------
% Initialise tissue rescaling weights
%__________________________________________________________________________
% Copyright (C) 2018 Wellcome Centre for Human Neuroimaging

    % =====================================================================
    % Read input from disk (if needed)
    [dat, datpath, opt] = fileToStruct(dat, opt);

    % =====================================================================
    % If the velocity is an observed -> nothing to do
    if defval(dat.v, '.observed', false)
        dat = structToFile(dat, datpath);
        return
    end

    % =====================================================================
    % Initial values
    dat.tpl.scale = ones(1,opt.model.nc);
              
    % =====================================================================
    % Exit
    dat = structToFile(dat, datpath);
end