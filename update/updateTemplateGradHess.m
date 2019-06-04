function dat = updateTemplateGradHess(dat, ~, opt)
% FORMAT dat = updateTemplateGradHess(dat, model, opt)
% dat   - Subject-specific data
% model - Model-specific data
% opt   - Options
%
% All inputs are structures that can either be in memory or on disk in the
% form of a mat file. In the latter case, it is read and, if needed,
% written back.
%--------------------------------------------------------------------------
% Update template gradient and hessian
%__________________________________________________________________________
% Copyright (C) 2018 Wellcome Centre for Human Neuroimaging

    % =====================================================================
    % Read input from disk (if needed)
    [dat, datpath, opt] = fileToStruct(dat, opt);
    
    % =====================================================================
    % Setup some constants
    if isfield(dat, 'model'),    noisemodel = dat.model;
    else,                        noisemodel = opt.model;  end
    
    % =====================================================================
    % Compute grad/hess
    [dat.tpl.g, dat.tpl.h] = ghTemplate(...
        noisemodel, ...            % Matching term (Categorical/Normal/...)
        dat.tpl.wmu, ...           % Warped+Softmaxed+Rescaled template
        dat.f.f, ...               % Observed image (class responsibilities)
        'ipsi',   dat.v.ipsi, ...  % Rigid+diffeo transform to push grad/hess
        'circ',  ~opt.tpl.bnd, ... % Boundary conditions 
        'lat',    opt.tpl.lat, ... % Template lattice (to push)
        'rotate', true, ...        % Rotate out null space
        'output', {dat.tpl.g, dat.tpl.h}, ... % Output file_array (saves memory)
        'par',    opt.par.within_subject, ... % parallelise processing? (usually no)
        'debug',  opt.ui.debug);   % Write debugging stuff? (usually no)
                             
    % =====================================================================
    % Exit
    dat = structToFile(dat, datpath);
end