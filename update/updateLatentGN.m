function dat = updateLatentGN(dat, model, opt)
% FORMAT dat = updateLatentGN(dat, model, opt)
% dat   - Subject-specific data
% model - Model-specific data
% opt   - Options
%
% All inputs are structures that can either be in memory or on disk in the
% form of a mat file. In the latter case, it is read and, if needed,
% written back.
%--------------------------------------------------------------------------
% Update latent parameters by Gauss-Newton
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
    if opt.iter.pena && dat.z.ok < 0
        dat.z.ok = dat.z.ok + 1;
        dat = structToFile(dat, datpath);
        return
    end
    
    % =====================================================================
    % Set a few constants
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
        [g, h] = ghMatchingLatent(...
            noisemodel, ...                     % Matching model (categorical/normal/...)
            dat.tpl.wmu, ...                    % Warped (+ softmaxed) template
            dat.f.f, ...                        % Observed matched image (responsibility)
            model.tpl.gmu, ...                  % (Log)-template spatial gradients
            model.pg.w, ...                     % Principal subspace
            'ipsi', dat.v.ipsi, ...             % Complete (rigid+diffeo) inverse transform
            'par', opt.par.within_subject, ...  % Parallelise stuff? (usually no)
            'debug', opt.ui.debug);             % Write debuging stuff? (usually no)

        % -----------------------------------------------------------------
        % Gradient/Hessian of the prior term
        [gz, hz] = ghPriorLatent(...
            dat.z.z, ...
            model.z.A + model.mixreg.w(2) * model.pg.ww, ...
            'debug', opt.ui.debug);
        g = g + gz;
        h = h + hz;
        clear gz hz

        % Part of shape-agnostic prior
        if model.mixreg.w(1) > 0 && checkarray(dat.v.r)
            m = spm_diffeo('vel2mom', single(numeric(dat.v.r)), double([opt.tpl.vs opt.pg.prm]));
            for k=1:opt.pg.K
                w1 = single(model.pg.w(:,:,:,:,k));
                g(k) = g(k) + model.mixreg.w(1) * w1(:)' * m(:);
            end
            clear w1
        end
        clear m
            
        % -----------------------------------------------------------------
        % Additional regularisation for robustness)
        h = spm_matcomp('LoadDiag', h);

        % -----------------------------------------------------------------
        % Compute search direction
        dz = -h\g;

        % -----------------------------------------------------------------
        % Line search
        result = lsLatent(...
            noisemodel, ...                     % Matching model (categorical/normal/...)
            dz, ....                            % Search direction
            dat.z.z, ...                        % Previous parameters
            dat.v.v, ...                        % Previous velocity
            dat.f.lb.val, ...                   % Previous matching term value
            model.pg.w, ...                     % Previous subspace
            model.tpl.a, ...                    % (Log)-template parameters
            dat.f.f, ...                        % Observed matched image (responsibility)
            'regz',    model.z.A, ...           % Prior precision matrix
            'geod',    model.mixreg.w(2), ...   % Weight of shape-agnostic prior
            'A',       A, ...                   % Affine transform
            'Mf',      dat.f.M, ...             % Image voxel-to-world
            'Mmu',     model.tpl.M, ...         % Template voxel-to-world 
            'nit',     lsiter, ...              % Number of line search iterations
            'itgr',    opt.iter.itg, ...        % Number of integration steps
            'prm',     opt.pg.prm, ...          % Differential operator parameters
            'bnd',     opt.pg.bnd, ...          % boundary conditions
            'par',     opt.par.within_subject, ... % Parallelise processing? (usually no)
            'verbose', opt.ui.verbose > 1, ...  % Talk during line search?
            'debug',   opt.ui.debug, ...        % Write debugging talk? (usually no)
            'pf',      dat.f.pf, ...            % File array to store the new pushed image
            'c',       dat.f.c, ...             % File array to store the new count image
            'wa',      dat.tpl.wa, ...          % File array to store the new warped log-template
            'wmu',     dat.tpl.wmu);            % File array to store the new warped+softmaxed template
            
        % -----------------------------------------------------------------
        % Store better values
        cumok = cumok || result.ok;
        compute_hessian = result.ok;
        if result.ok
            dat.z.z       = result.z;
            dat.z.zz      = result.z * result.z';
            dat.f.lb.val  = result.llm;
            dat.v.ipsi    = copyarray(result.ipsi, dat.v.ipsi);
            dat.v.v       = copyarray(result.v,    dat.v.v);
            if model.mixreg.w(2)
                m = spm_diffeo('vel2mom', result.v, double([opt.tpl.vs opt.pg.prm]));
                dat.v.lb.regv = result.v(:)' * m(:);
                clear m
            end
            if strcmpi(opt.tpl.update, 'ml')
                dat.f.pf      = copyarray(result.pf,   dat.f.pf);
                dat.f.c       = copyarray(result.c,    dat.f.c);
            else
                rmarray(result.pf);
                rmarray(result.c);
            end
            dat.f.bb      = result.bb;
            dat.tpl.wmu   = copyarray(result.wmu, dat.tpl.wmu);
            rmarray(result.wa);
        else
            break
        end

    end % < GN iterations
    
    % =====================================================================
    % Don't try next time if it failed
    if cumok
        dat.z.ok2 = 0;
        dat.z.ok  = 1; 
    else
        if opt.iter.pena
            dat.z.ok2 = dat.z.ok2 - 1;
            dat.z.ok  = dat.z.ok2; 
        else
            dat.z.ok2 = 0;
            dat.z.ok  = 0;
        end
    end
    
    % =====================================================================
    % Update Hessian for Laplace approximation
    if compute_hessian

        % Likelihood part
        h = ghMatchingLatent(...
            noisemodel, ...                     % Matching model (categorical/normal/...)
            dat.tpl.wmu, ...                    % Warped (+ softmaxed) template
            dat.f.f, ...                        % Observed matched image (responsibility)
            model.tpl.gmu, ...                  % (Log)-template spatial gradients
            model.pg.w, ...                     % Principal subspace
            'ipsi',    dat.v.ipsi, ...          % Complete (rigid+diffeo) inverse transform
            'hessian', true, ...                % Compute only hessia
            'par',     opt.par.within_subject, ... % Parallelise stuff? (usually no)
            'debug',   opt.ui.debug);           % Write debuging stuff? (usually no)

        % prior part
        [~, hz] = ghPriorLatent(...
            dat.z.z, ...
            model.z.A + model.mixreg.w(2) * model.pg.ww, ...
            'debug', opt.ui.debug);
        h = h + hz;
        clear hz

        % Additional regularisation for robustness
        h = spm_matcomp('LoadDiag', h);
    end
    dat.z.S = spm_matcomp('Inv', h);
              
    % =====================================================================
    % Exit
    dat = structToFile(dat, datpath);
end