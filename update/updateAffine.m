function dat = updateAffine(dat, model, opt)
% FORMAT dat = updateAffine(dat, model, opt)
% dat   - Subject-specific data
% model - Model-specific data
% opt   - Options
%
% All inputs are structures that can either be in memory or on disk in the
% form of a mat file. In the latter case, it is read and, if needed,
% written back.
%--------------------------------------------------------------------------
% Update affine parameters by Gauss-Newton
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
    if opt.iter.pena && dat.q.ok < 0
        dat.q.ok = dat.q.ok + 1;
        dat = structToFile(dat, datpath);
        return
    end

    % =====================================================================
    % Load stuff if needed
    if isempty(defval(dat.buffer, 'f.f', []))
        if opt.buf, dat.buffer.f.f = numeric(dat.f.f);
        else,       dat.buffer.f.f = dat.f.f;                   end
    end
    if isempty(defval(dat.buffer, 'tpl.wmu', []))
        if opt.buf, dat.buffer.tpl.wmu = numeric(dat.tpl.wmu);
        else,       dat.buffer.tpl.wmu = dat.tpl.wmu;           end
    end
    if isempty(defval(dat.buffer, 'tpl.a', []))
        if opt.buf, dat.buffer.tpl.a = numeric(model.tpl.a);
        else,       dat.buffer.tpl.a = model.tpl.a;             end
    end
    if isempty(defval(dat.buffer, 'tpl.gmu', []))
        if opt.buf, dat.buffer.tpl.gmu = numeric(model.tpl.gmu);
        else,       dat.buffer.tpl.gmu = model.tpl.gmu;         end
    end
    if isempty(defval(dat.buffer, 'v.ipsi', []))
        if opt.buf, dat.buffer.v.ipsi = numeric(dat.v.ipsi);
        else,       dat.buffer.v.ipsi = dat.v.ipsi;             end
    end
    if isfield(dat.v, 'v') && isempty(defval(dat.buffer, 'v.v', []))
        if opt.buf, dat.buffer.v.v = numeric(dat.v.v);
        else,       dat.buffer.v.v = dat.v.v;                   end
    end
    
    % =====================================================================
    % If non-rigid transformation
    % > compute direct diffeomorphism and jacobian determinant
    if isfield(dat.v, 'v')
        [iphi, phi, jac] = exponentiateVelocity(...
            dat.buffer.v.v, ...        % Current velocity
            'iphi', 'phi', 'jac', ...  % What to compute?
            'itgr',  opt.iter.itg, ... % Nb of integration steps
            'vs',    opt.tpl.vs, ...   % Velocity voxel size
            'prm',   opt.pg.prm, ...   % Registration regularisation param
            'bnd',   opt.pg.bnd, ...   % Boundary condition
            'debug', opt.ui.debug);    % Write debugging stuff (usually no)
        if opt.model.dim == 2
            iphi(:,:,:,3) = 1;
            phi(:,:,:,3)  = 1;
        end
    else
        % Pure rigid/affine case
        iphi = spm_warps('identity', opt.tpl.lat);
        phi  = [];
        jac  = [];
    end
    
    % =====================================================================
    % Set a few constants
    if isfield(opt.q, 'gniter'), gniter = opt.q.gniter;
    else,                        gniter = opt.iter.gn;    end
    if isfield(opt.q, 'lsiter'), lsiter = opt.q.lsiter;
    else,                        lsiter = opt.iter.ls;    end
    if isfield(dat, 'model'),    noisemodel = dat.model;
    else,                        noisemodel = opt.model;  end
    
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
        [g, h] = ghMatchingAffine(...
            noisemodel, ...                  % Matching model (categorical/normal/...)
            dat.buffer.tpl.wmu, ...          % Warped (+ softmaxed) template
            dat.buffer.f.f, ...              % Observed matched image (responsibility)
            dat.buffer.tpl.gmu, ...          % (Log)-template spatial gradients
            dat.q.A, ...                     % Current rigid/affine transform
            opt.q.B, ...                     % Rigid/affine Lie basis
            phi, ...                         % Direct diffeomorphism
            jac, ...                         % Jacobian determinant of the direct diffeo
            'ipsi',   dat.buffer.v.ipsi, ... % Complete (rigid+diffeo) inverse transform
            'circ',  ~opt.tpl.bnd, ...       % Boundary conditions
            'Mmu',    model.tpl.M, ...       % Template voxel-to-world
            'approx', opt.q.hapx, ....       % Approximate hessian? (usually true)
            'par',    opt.par.within_subject, ... % Parallelise stuff? (usually no)
            'debug',  opt.ui.debug);          % Write debuging stuff? (usually no)

        % -----------------------------------------------------------------
        % Gradient/Hessian of the prior term (only if not rigid)
        if opt.q.Mr
            rind = opt.q.rind;
            [gq, hq] = ghPriorAffine(dat.q.q(rind), model.q.A, 'debug', opt.ui.debug);
            g(rind)      = g(rind)      + gq;
            h(rind,rind) = h(rind,rind) + hq;
            clear gq hq
            A = model.q.A;
        else
            rind = [];
            A    = [];
        end

        % -----------------------------------------------------------------
        % Additional regularisation for robustness)
        h = spm_matcomp('LoadDiag', h);

        % -----------------------------------------------------------------
        % Compute search direction
        dq = -h\g;

        % -----------------------------------------------------------------
        % Line search
        result = lsAffine(...
            noisemodel, ...             % Matching model (categorical/normal/...)
            dq, ...                     % Search direction
            dat.q.q, ...                % Previous parameters
            dat.f.lb.val, ...           % Previous matching term value
            dat.buffer.tpl.a, ...       % (Log)-template parameters
            dat.f.f, ...                % Observed matched image (responsibility)
            'B',       opt.q.B, ...     % Rigid/affine Lie basis
            'regq',    A, ...           % Prior precision matrix (only if not rigid)
            'rind',    rind, ...        % Indices of regularised (affine) param
            'iphi',    iphi, ...        % Inverse diffeomorphism
            'Mf',      dat.f.M, ...     % Image voxel-to-world
            'Mmu',     model.tpl.M, ... % Template voxel-to-world 
            'nit',     lsiter, ...      % Number of line search iterations
            'par',     opt.par.within_subject, ... % Parallelise processing? (usually no)
            'verbose', opt.ui.verbose > 1, ...     % Talk during line search?
            'debug',   opt.ui.debug);   % Write debugging talk? (usually no)
%             'pf',      dat.f.pf, ...    % File array to store the new pushed image
%             'c',       dat.f.c, ...     % File array to store the new count image
%             'wa',      dat.tpl.wa, ...  % File array to store the new warped log-template
%             'wmu',     dat.tpl.wmu);    % File array to store the new warped+softmaxed template

        % -----------------------------------------------------------------
        % Store better values
        cumok = cumok || result.ok;
        compute_hessian = result.ok;
        if result.ok
            dat.q.q       = result.q;
            dat.q.qq      = result.q * result.q';
            dat.q.A       = result.A;
            dat.f.lb.val  = result.llm;
            dat.v.ipsi    = copyarray(result.ipsi, dat.v.ipsi);
            if opt.buf, dat.buffer.v.ipsi = result.ipsi; end
            if strcmpi(opt.tpl.update, 'ml')
                dat.f.pf      = copyarray(result.pf,   dat.f.pf);
                dat.f.c       = copyarray(result.c,    dat.f.c);
            else
                rmarray(result.pf);
                rmarray(result.c);
            end
            dat.f.bb      = result.bb;
            dat.tpl.wmu   = copyarray(result.wmu, dat.tpl.wmu);
            if opt.buf, dat.buffer.tpl.wmu = result.wmu; end
            rmarray(result.wa);
        else
            break
        end

    end % < GN iterations
    
    % =====================================================================
    % Don't try next time if it failed
    if cumok
        dat.q.ok2 = 0;
        dat.q.ok  = 1; 
    else
        if opt.iter.pena
            dat.q.ok2 = dat.q.ok2 - 1;
            dat.q.ok  = dat.q.ok2; 
        else
            dat.q.ok2 = 0;
            dat.q.ok  = 0;
        end
    end
    
    % =====================================================================
    % Update Hessian for Laplace approximation
    if opt.q.Mr
        if compute_hessian

            % Likelihood part
            h = ghMatchingAffine(...
                noisemodel, ...             % Matching model (categorical/normal/...)
                dat.buffer.tpl.wmu, ...     % Warped (+ softmaxed) template
                dat.buffer.f.f, ...         % Observed matched image (responsibility)
                dat.buffer.tpl.gmu, ...     % (Log)-template spatial gradients
                dat.q.A, ...                % Current rigid/affine transform
                opt.q.B, ...                % Rigid/affine Lie basis
                phi, ...                    % Direct diffeomorphism
                jac, ...                    % Jacobian determinant of the direct diffeo
                'ipsi',    dat.buffer.v.ipsi, ... % Complete (rigid+diffeo) inverse transform
                'circ',  ~opt.tpl.bnd, ...  % Boundary conditions
                'hessian', true, ...        % Do not compute gradient
                'Mmu',     model.tpl.M, ... % Template voxel-to-world
                'approx',  opt.q.hapx, .... % Approximate hessian? (usually true)
                'par',     opt.par.within_subject, ... % Parallelise stuff? (usually no)
                'debug',   opt.ui.debug);   % Write debuging stuff? (usually no)

            % prior part
            [~, hq] = ghPriorAffine(dat.q(rind), model.q.A, 'debug', opt.ui.debug);
            h(rind,rind) = h(rind,rind) + hq;
            clear hq

            % Additional regularisation for robustness
            h = spm_matcomp('LoadDiag', h);
        end
        dat.q.S = spm_matcomp('Inv', h);
    else
        dat.q.S = 0;
    end
              
    % =====================================================================
    % Exit
    dat = structToFile(dat, datpath);
end