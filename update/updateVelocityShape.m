function dat = updateVelocityShape(dat, model, opt)
% FORMAT dat = updateVelocityShape(dat, model, opt)
% dat   - Subject-specific data
% model - Model-specific data
% opt   - Options
%
% All inputs are structures that can either be in memory or on disk in the
% form of a mat file. In the latter case, it is read and, if needed,
% written back.
%--------------------------------------------------------------------------
% Update velocity field (shape model version)
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
    % Penalise previous failure
    % > if option activated and previous update failed, do not try
    if opt.iter.pena && dat.v.ok < 0
        dat.v.ok = dat.v.ok + 1;
        dat = structToFile(dat, datpath);
        return
    end
    
    % =====================================================================
    % Set a few constants
    spm_diffeo('boundary', opt.pg.bnd);
    if isfield(opt.q, 'gniter'), gniter = opt.q.gniter;
    else,                        gniter = opt.iter.gn;    end
    if isfield(opt.q, 'lsiter'), lsiter = opt.q.lsiter;
    else,                        lsiter = opt.iter.ls;    end
    if isfield(dat, 'model'),    noisemodel = dat.model;
    else,                        noisemodel = opt.model;  end
    if isfield(dat, 'q') && isfield(dat.q, 'A'),  A = dat.q.A;
    else,                                         A = eye(4);  end
    
    % =====================================================================
    % Gauss-Newton iterations
    % It is useful to actually find a mode of the posterior (and not 
    % only an improved value) when we use the Laplace precision for  
    % the update of W. In that case, setting gnit > 1 might help  
    % converge faster.
    cumok = false;
    for i=1:gniter

        % -----------------------------------------------------------------
        % Gradient/Hessian of the likelihood term
        [g, h] = ghMatchingVel(...
            noisemodel, ...             % Matching model (categorical/normal/...)
            dat.tpl.wmu, ...            % Warped (+ softmaxed) template
            dat.f.f, ...                % Observed matched image (responsibility)
            model.tpl.gmu, ...          % (Log)-template spatial gradients
            'ipsi',  dat.v.ipsi, ...    % Complete (rigid+diffeo) inverse transform
            'circ',  ~opt.tpl.bnd, ...  % Boundary conditions
            'par',   opt.par.within_subject, ... % Parallelise stuff? (usually no)
            'debug', opt.ui.debug);     % Write debuging stuff? (usually no)


        % -----------------------------------------------------------------
        % Gradient/Hessian of the prior term
        % g = g - w * lam * LWz
        wz = reconstructVelocity(...
                'latent',   dat.z.z, ...
                'subspace', model.pg.w, ...
                'par',      opt.par.within_subject);
        g = g - spm_diffeo('vel2mom', wz, double([opt.tpl.vs (model.mixreg.w(1) * model.v.l * opt.pg.prm)]));
        v = numeric(dat.v.v);
        r = v - wz;
        clear wz
        % g = g + (w * lam + (1-w)) * Lv
        g = g + ghPriorVel(v, opt.tpl.vs, (model.mixreg.w(1)*model.v.l + model.mixreg.w(2)) * opt.pg.prm, opt.pg.bnd);
        if opt.model.dim == 2
            g(:,:,:,3) = 0;   % 2D case: ensure null gradient in 3rd dim
        end

        % -----------------------------------------------------------------
        % Search direction
        dv = -spm_diffeo('fmg', single(h), single(g), ...
            double([opt.tpl.vs (model.mixreg.w(1)*model.v.l + model.mixreg.w(2)) * opt.pg.prm 2 2]));
        clear g
        if opt.model.dim == 2
            dv(:,:,:,3) = 0;   % 2D case: ensure null velocity in 3rd dim
        end

        % -----------------------------------------------------------------
        % Line search
        result = lsVelocityShape(...
            noisemodel, ...             % Matching model (categorical/normal/...)
            dv, ...                     % Search direction
            r, ...                      % Previous residual field
            dat.f.lb.val, ...
            model.tpl.a, ...
            dat.f.f, ...
            'v0',       v, ...                 % Previous velocity field (only needed if w(2) > 0)
            'lam',      model.v.l, ...         % Residual precision
            'w',        model.mixreg.w(2), ... % Mixture of regularisations weight
            'prm',      opt.pg.prm, ...        % Regularisation parameters
            'itgr',     opt.iter.itg, ...      % Number of integration steps
            'bnd',      opt.pg.bnd, ...        % Boundary conditions
            'A',        A, ...                 % Rigid/affine transform
            'Mf',       dat.f.M, ...           % Image voxel-to-world
            'Mmu',      model.tpl.M, ...       % Template voxel-to-world
            'nit',      lsiter,  ...           % Line search iterations
            'par',      opt.par.within_subject, ... % Parallelise processing? (usually no)
            'verbose',  opt.ui.verbose>1, ...  % Talk during line search?
            'debug',    opt.ui.debug, ...      % Write debugging talk? (usually no)
            'pf',       dat.f.pf, ...          % File array to store the new pushed image
            'c',        dat.f.c, ...           % File array to store the new count image
            'wa',       dat.tpl.wa, ...        % File array to store the new warped log-template
            'wmu',      dat.tpl.wmu);          % File array to store the new warped+softmaxed template

        % -----------------------------------------------------------------
        % Store better values
        cumok = cumok || result.ok;
        compute_hessian = result.ok;
        if result.ok
            dat.f.lb.val  = result.match;
            dat.v.v       = copyarray(result.v,    dat.v.v);
            m = spm_diffeo('vel2mom', result.v, double([opt.tpl.vs opt.pg.prm]));
            if model.mixreg.w(2)
                dat.v.lb.regv = result.v(:)' * m(:);
            end
            clear m
            dat.v.r       = copyarray(result.r,    dat.v.r);
            dat.v.ipsi    = copyarray(result.ipsi, dat.v.ipsi);
            dat.f.pf      = copyarray(result.pf,   dat.f.pf);
            dat.f.c       = copyarray(result.c,    dat.f.c);
            dat.f.bb      = result.bb;
            dat.tpl.wmu   = copyarray(result.wmu, dat.tpl.wmu);
            rmarray(result.wa);
            r = result.r;
            ipsi = dat.v.ipsi;
        else
            break
        end
        clear result

    end % < GN iterations
    
    
    % =====================================================================
    % Don't try next time if it failed
    if cumok
        dat.v.ok2 = 0;
        dat.v.ok  = 1; 
    else
        if opt.iter.pena
            dat.v.ok2 = dat.v.ok2 - 1;
            dat.v.ok  = dat.v.ok2; 
        else
            dat.v.ok2 = 0;
            dat.v.ok  = 0;
        end
    end

    % =====================================================================
    % Lower bound elements
    
    % ---------------------------------------------------------------------
    % Residual regularisation
    m = spm_diffeo('vel2mom', single(numeric(r)), double([opt.tpl.vs opt.pg.prm]));
    dat.v.lb.reg = r(:)' * m(:);
    clear r m

    % ---------------------------------------------------------------------
    % Update Hessian for Laplace approximation
    if compute_hessian
        h = ghMatchingVel(...
            noisemodel, ...             % Matching model (categorical/normal/...)
            dat.tpl.wmu, ...            % Warped (+ softmaxed) template
            dat.f.f, ...                % Observed matched image (responsibility)
            model.tpl.gmu, ...          % (Log)-template spatial gradients
            'ipsi',    ipsi, ...        % Complete (rigid+diffeo) inverse transform
            'circ',  ~opt.tpl.bnd, ...  % Boundary conditions
            'hessian', true, ...        % Do not compute gradient
            'par',     opt.par.within_subject, ... % Parallelise stuff? (usually no)
            'debug',   opt.ui.debug);   % Write debuging stuff? (usually no)
    end
    
    % ---------------------------------------------------------------------
    % Tr((H+L)\L)
    dat.v.lb.tr = trapprox((model.mixreg.w(1)*model.v.l + model.mixreg.w(2))*opt.pg.prm, h, 'vs', opt.tpl.vs);
    dat.v.lb.tr = dat.v.lb.tr / (model.mixreg.w(1)*model.v.l + model.mixreg.w(2));

    % ---------------------------------------------------------------------
    % LogDet(H+L)
    dat.v.lb.ld = ldapprox((model.mixreg.w(1)*model.v.l + model.mixreg.w(2))*opt.pg.prm, h, 'vs', opt.tpl.vs);
    clear h
                      
    % =====================================================================
    % Exit
    dat = structToFile(dat, datpath);
end