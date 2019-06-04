function dat = initVelocityShapeLaplace(dat, model, opt)
% FORMAT dat = initVelocityShapeLaplace(dat, model, opt)
% dat   - Subject-specific data
% model - Model-specific data
% opt   - Options
%
% All inputs are structures that can either be in memory or on disk in the
% form of a mat file. In the latter case, it is read and, if needed,
% written back.
%--------------------------------------------------------------------------
% Initialise velocity Laplace approximation (shape model version)
%__________________________________________________________________________
% Copyright (C) 2018 Wellcome Centre for Human Neuroimaging

    % =====================================================================
    % Read input from disk (if needed)
    [dat, datpath, model, ~, opt] = fileToStruct(dat, model, opt);
    
    % =====================================================================
    % If the velocity is an observed -> nothing to do
    if defval(dat.v, '.observed', false)
        dat = structToFile(dat, datpath);
        return
    end

    % =====================================================================
    % Lower bound elements
    spm_diffeo('boundary', opt.pg.bnd);
    
    % ---------------------------------------------------------------------
    % Residual regularisation
    dat.v.lb.reg  = 0;
    dat.v.lb.regv = 0;

    % ---------------------------------------------------------------------
    % Update Hessian for Laplace approximation
    if isfield(dat, 'model'),    noisemodel = dat.model;
    else,                        noisemodel = opt.model;  end
    h = ghMatchingVel(...
        noisemodel, ...             % Matching model (categorical/normal/...)
        dat.tpl.wmu, ...            % Warped (+ softmaxed) template
        dat.f.f, ...                % Observed matched image (responsibility)
        model.tpl.gmu, ...          % (Log)-template spatial gradients
        'ipsi',    dat.v.ipsi, ...  % Complete (rigid+diffeo) inverse transform
        'hessian', true, ...        % Do not compute gradient
        'par',     opt.par.within_subject, ... % Parallelise stuff? (usually no)
        'debug',   opt.ui.debug);   % Write debuging stuff? (usually no)
    
    % ---------------------------------------------------------------------
    % Tr((H+L)\L)
    dat.v.lb.tr = trapprox((model.mixreg.w(1)*model.v.l + model.mixreg.w(2))*opt.pg.prm, h, 'vs', opt.tpl.vs);
    dat.v.lb.tr = dat.v.lb.tr / (model.mixreg.w(1)*model.v.l + model.mixreg.w(2));

    % ---------------------------------------------------------------------
    % LogDet(H+L)
    dat.v.lb.ld = ldapprox((model.mixreg.w(1)*model.v.l + model.mixreg.w(2))*opt.pg.prm, h, 'vs', opt.tpl.vs);
    clear h
    
    % ---------------------------------------------------------------------
    dat.v.lb.uncty = trace(dat.z.S*model.pg.ww);
    
    % =====================================================================
    % Exit
    dat = structToFile(dat, datpath);
end