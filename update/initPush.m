function dat = initPush(dat, ~, opt)
% FORMAT dat = initPush(dat, model, opt)
% dat   - Subject-specific data
% model - Model-specific data
% opt   - Options
%
% All inputs are structures that can either be in memory or on disk in the
% form of a mat file. In the latter case, it is read and, if needed,
% written back.
%--------------------------------------------------------------------------
% Initialise pushed images
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
    [dat.f.pf, dat.f.c, dat.f.bb] = pushImage(...
        dat.v.ipsi, ...                       % Complete transform
        dat.f.f, ...                          % Original image
        opt.tpl.lat, ...                      % Template lattice
        'output', {dat.f.pf, dat.f.c}, ...    % Output file array (save memory)
        'par',    opt.par.within_subject, ... % Parallelise processing? (usually no)
        'debug',  opt.ui.debug);              % Write debugging stuff? (usually no)
              
    % =====================================================================
    % Exit
    dat = structToFile(dat, datpath);
end