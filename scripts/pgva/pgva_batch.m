function varargout = pgva_batch(id, varargin)
% _________________________________________________________________________
%
%             Batch functions for the PGVA model (2018)
% _________________________________________________________________________
%
% Collection of functions that require updating some subject-specific
% variables. Depending on the needs, as much as possible is done in order
% to minimise i/o, memory and computational costs.
%
% --------------
% Initialisation
% --------------
%
% FORMAT dat = pgva_batch('InitPush', dat, model, opt)
%   Initialise warped images (Change of lattice, but no deformation)
%
% FORMAT [dat, model] = pgva_batch('InitPull', dat, model, opt)
%   Initialise warped template (Change of lattice, but no deformation)
%   + matching log-likelihood
%
% FORMAT [dat, model] = pgva_batch('InitAffine', mode, dat, model, opt)
%   Initialise affine coordinates ('zero'/'rand') + lower bound
%
% FORMAT [dat, model] = pgva_batch('InitVelocity', mode, dat, model, opt)
%   Initialise velocity fields ('zero'/'rand') + lower bound
%
% FORMAT [dat, model] = pgva_batch('InitLatent', mode, dat, model, opt)
%   Initialise latent coordinates ('zero'/'rand') + lower bound
%   If random, insure they are zero-centered and orthogonal.
%
% ---
% Fit
% ---
%
% FORMAT [dat, model] = pgva_batch('FitAffine', dat, model, opt)
%   Gauss-Newton update of the affine coordinates
%
% FORMAT [dat, model] = pgva_batch('FitLatent', dat, model, opt)
%   Gauss-Newton update of the latent coordinates
%
% FORMAT [dat, model] = pgva_batch('FitVelocity', dat, model, opt)
%   Gauss-Newton update of the residual fields
%
% ----------
% Correction
% ----------
% 
% FORMAT dat = pgva_batch('RotateSubspace', R, iR, dat, model, opt)
%   Apply a rotation matrix to all latent coordinates
%
% % FORMAT dat = pgva_batch('CentreLatent', dat, opt, (mean))
% %   Subtract the mean to all latent coordinates
%
% ------
% Update
% ------
%
% FORMAT [dat, model] = pgva_batch('LB', var, dat, model, opt)
%   Update or compute the part of the lower bound that depends on var.
%   var can be 'Matching'/'Lambda'/'Subspace'.
%
% FORMAT dat = pgva_batch('Momentum', dat, opt)
%   Compute momentum from observed velocities
% _________________________________________________________________________
%
% The following subfunctions should not be used usually. They are only
% needed when distributing jobs on a cluster.
% 
% FORMAT dat = pgva_batch('OneInitPush',     dat, model, opt)
% FORMAT dat = pgva_batch('OneInitPull',     dat, model, opt)
% FORMAT dat = pgva_batch('OneInitVelocity', dat, model, opt, mode)
% FORMAT dat = pgva_batch('OneFitAffine',    dat, model, opt)
% FORMAT dat = pgva_batch('OneFitLatent',    dat, model, opt)
% FORMAT dat = pgva_batch('OneFitVelocity',  dat, model, opt)
% FORMAT dat = pgva_batch('OneLB',           dat, model, opt, var)
% _________________________________________________________________________

    switch lower(id)
        % -------
        %  BATCH
        % -------
        case 'initpush'
            [varargout{1:nargout}] = batchInitPush(varargin{:});
        case 'initpull'
            [varargout{1:nargout}] = batchInitPull(varargin{:});
        case 'initaffine'
            [varargout{1:nargout}] = batchInitAffine(varargin{:});
        case 'initlatent'
            [varargout{1:nargout}] = batchInitLatent(varargin{:});
        case 'initvelocity'
            [varargout{1:nargout}] = batchInitVelocity(varargin{:});
        case 'rotatesubspace'
            [varargout{1:nargout}] = batchRotateSubspace(varargin{:});
        case 'fitaffine'
            [varargout{1:nargout}] = batchFitAffine(varargin{:});
        case 'fitlatent'
            [varargout{1:nargout}] = batchFitLatent(varargin{:});
        case 'fitvelocity'
            [varargout{1:nargout}] = batchFitVelocity(varargin{:});
        case 'lb'
            [varargout{1:nargout}] = batchLB(varargin{:});
        case 'momentum'
            [varargout{1:nargout}] = batchMomentum(varargin{:});
        % ------
        %  STEP
        % ------
        case 'oneinitpush'
            [varargout{1:nargout}] = oneInitPush(varargin{:});
        case 'oneinitpull'
            [varargout{1:nargout}] = oneInitPull(varargin{:});
        case 'oneinitvelocity'
            [varargout{1:nargout}] = oneInitVelocity(varargin{:});
        case 'onefitaffine'
            [varargout{1:nargout}] = oneFitAffine(varargin{:});
        case 'onefitlatent'
            [varargout{1:nargout}] = oneFitLatent(varargin{:});
        case 'onefitvelocity'
            [varargout{1:nargout}] = oneFitVelocity(varargin{:});
        case 'onelb'
            [varargout{1:nargout}] = oneLB(varargin{:});
        case 'onemomentum'
            [varargout{1:nargout}] = oneMomentum(varargin{:});
    end

end

% -------------------------------------------------------------------------
%    Helper
% -------------------------------------------------------------------------

function before = plotBatchBegin(name)
    before = 0;
    fprintf('%10s | %10s | ', 'Batch', name);
    tic;
end

function before = plotBatch(i, batch, n, total, before)
    exact = min(n, i*batch)/n*total;
    trunc = floor(exact);
    step  = trunc - before;
    before = trunc;
    fprintf(repmat('.', 1, step));
end


function plotBatchEnd
    fprintf(' | %6.3fs\n', toc);
end

% -------------------------------------------------------------------------
%    Momentum
% -------------------------------------------------------------------------

function dat = oneMomentum(dat, ~, opt)

    if ~dat.v.observed
        return
    end
    
    m = spm_diffeo('vel2mom', single(numeric(dat.v.v)), double([opt.tpl.vs opt.pg.prm]));
    dat.v.m = prepareOnDisk(dat.v.m, size(m));
    dat.v.m(:) = m(:);
    
end

function dat = batchMomentum(dat, opt)

    % --- Detect parallelisation scheme
    if strcmpi(opt.split.loop, 'subject') && opt.split.par > 0
        batch = opt.split.batch;
    else
        batch = 1;
    end
    
    % --- Batch processing
    N = numel(dat);
    if opt.ui.verbose, before = plotBatchBegin('Init Mom'); end
    for i=1:ceil(N/batch)
        n1 = (i-1)*batch + 1;
        ne = min(N, i*batch);

        if opt.ui.verbose, before = plotBatch(i, batch, N, 50, before); end
    
        dat(n1:ne) = distribute('OneMomentum', dat(n1:ne), struct, opt);
        
    end
    if opt.ui.verbose, plotBatchEnd; end
end

% -------------------------------------------------------------------------
%    Init Push
% -------------------------------------------------------------------------

function dat = oneInitPush(dat, model, opt)

    if ~dat.f.observed
        return
    end

    % --- Detect parallelisation scheme
    if strcmpi(opt.split.loop, 'subject')
        loop = '';
        par  = 0;
    else
        loop = opt.split.loop;
        par  = opt.split.par;
    end
    
    % Warp
    % ----
    iphi = spm_warps('identity', opt.tpl.lat);
    ipsi = reconstructIPsi(eye(4), iphi, 'lat', opt.tpl.lat, ...
                           'Mf',  dat.f.M, 'Mmu', model.tpl.M, ...
                           'debug', opt.ui.debug);
    clear iphi
    [dat.f.pf, dat.f.c, dat.f.bb] = pushImage(ipsi, dat.f.f, opt.tpl.lat, ...
                                              'output', {dat.f.pf, dat.f.c}, ...
                                              'loop', loop, 'par', par, ...
                                              'debug', opt.ui.debug);
    clear ipsi

end

function dat = batchInitPush(dat, model, opt)

    % --- Detect parallelisation scheme
    if strcmpi(opt.split.loop, 'subject') && opt.split.par > 0
        batch = opt.split.batch;
    else
        batch = 1;
    end
    
    % --- Batch processing
    N = numel(dat);
    if opt.ui.verbose, before = plotBatchBegin('Init Push'); end
    for i=1:ceil(N/batch)
        n1 = (i-1)*batch + 1;
        ne = min(N, i*batch);

        if opt.ui.verbose, before = plotBatch(i, batch, N, 50, before); end
    
        dat(n1:ne) = distribute('OneInitPush', dat(n1:ne), model, opt);
        
    end
    if opt.ui.verbose, plotBatchEnd; end
end

% -------------------------------------------------------------------------
%    Init Pull
% -------------------------------------------------------------------------

function dat = oneInitPull(dat, model, opt)

    if ~dat.f.observed
        return
    end

    % --- Detect parallelisation scheme
    if strcmpi(opt.split.loop, 'subject')
        loop = '';
        par  = 0;
    else
        loop = opt.split.loop;
        par  = opt.split.par;
    end
    
    % --- Detect noise model
    if isfield(dat, 'model')
        noisemodel = dat.model;
    else
        noisemodel = opt.model;
    end
    
    % Warp
    % ----
    iphi = spm_warps('identity', opt.tpl.lat);
    ipsi = reconstructIPsi(eye(4), iphi, 'lat', opt.tpl.lat, ...
                           'Mf',  dat.f.M, 'Mmu', model.tpl.M, ...
                           'debug', opt.ui.debug);
    clear iphi
    if opt.tpl.cat
        dat.tpl.wa = warp(ipsi, model.tpl.a, opt.tpl.itrp, opt.tpl.bnd, ...
                          'par', par, 'output', dat.tpl.wa);
        dat.tpl.wmu = reconstructProbaTemplate(dat.tpl.wa, ...
                                               'output', dat.tpl.wmu, ...
                                               'loop', loop, 'par', par, ...
                                               'debug', opt.ui.debug);
        dat.tpl.wa = rmarray(dat.tpl.wa);
    else
        dat.tpl.wmu = warp(ipsi, model.tpl.mu, opt.tpl.itrp, opt.tpl.bnd, ...
                           'par', par, 'output', dat.tpl.wmu, ...
                           'debug', opt.ui.debug);
    end
    clear ipsi
    
    % Matching likelihood
    % -------------------
    dat.f.lb.type = 'll';
%     dat.f.lb.val  = llMatching(noisemodel, dat.tpl.wmu, dat.f.f, ...
%                                'loop', loop, 'par', par);
    dat.f.lb.val  = llMatching(noisemodel, model.tpl.mu, dat.f.pf, dat.f.c, ...
                               'loop', loop, 'par', par, 'bb', dat.f.bb);

end

function [dat, model] = batchInitPull(dat, model, opt)

    % --- Detect parallelisation scheme
    if strcmpi(opt.split.loop, 'subject') && opt.split.par > 0
        batch = opt.split.batch;
    else
        batch = 1;
    end
    
    % --- Init lower bound
    model.lb.m.val = 0;
    
    % --- Batch processing
    N = numel(dat);
    if opt.ui.verbose, before = plotBatchBegin('Init Pull'); end
    for i=1:ceil(N/batch)
        n1 = (i-1)*batch + 1;
        ne = min(N, i*batch);

        if opt.ui.verbose, before = plotBatch(i, batch, N, 50, before); end
    
        dat(n1:ne) = distribute('OneInitPull', dat(n1:ne), model, opt);
        
        for n=n1:ne
            if dat(n).f.observed
                model.lb.m.val = model.lb.m.val + dat(n).f.lb.val;
            end
        end
    end
    if opt.ui.verbose, plotBatchEnd; end
end

% -------------------------------------------------------------------------
%    Init Affine
% -------------------------------------------------------------------------

function [dat, model] = batchInitAffine(mode, dat, model, opt)

    % --- Don't parallelise here, there is no need
    
    switch lower(mode)
        case 'zero'
            create = @zeros;
        case 'rand'
            create = @randn;
        otherwise
            error('Unknown mode %s', mode)
    end
    
    % Init model
    % ----------
    model.q.q  = zeros(opt.q.M, 1);
    model.q.qq = zeros(opt.q.M);
    model.q.S  = zeros(opt.q.M);
    model.lb.q.val = 0;
    
    % Init subjects
    % -------------
    nq = 0;
    if opt.ui.verbose, before = plotBatchBegin('Init Q'); end
    for n=1:numel(dat)
        if dat(n).f.observed
            nq = nq + 1;
            if opt.ui.verbose, before = plotBatch(nq, 1, opt.f.N, 50, before); end
            q           = create([opt.q.M, 1]);
            dat(n).q.q  = q;
            dat(n).q.qq = q*q';
            dat(n).q.S  = zeros(opt.q.M);
            dat(n).q.A  = exponentiateAffine(q, opt.q.B);
            model.q.q   = model.q.q  + dat(n).q.q;
            model.q.qq  = model.q.qq + dat(n).q.qq;
            dat(n).q.ok  = 1; % for GN failure tracking
            dat(n).q.ok2 = 0; % for GN failure tracking
            if opt.q.Mr
                rind = opt.q.rind;
                dat(n).q.lb.type = 'kl';
                dat(n).q.lb.val  = -0.5*( trace((dat(n).q.S(rind,rind) + dat(n).q.qq(rind,rind)) * model.q.A) ...
                                          - spm_matcomp('LogDet', model.q.A) ...
                                          - spm_matcomp('LogDet', dat(n).q.S(rind,rind)) ...
                                          - opt.q.Mr );
            else
                dat(n).q.lb.type = '';
                dat(n).q.lb.val = 0;
            end
            model.lb.q.val = model.lb.q.val + dat(n).q.lb.val;
        end
    end
    if opt.ui.verbose, plotBatchEnd; end

end

% -------------------------------------------------------------------------
%    Init Velocity
% -------------------------------------------------------------------------

function dat = oneInitVelocity(dat, model, opt, mode)

    % --- Detect parallelisation scheme
    if strcmpi(opt.split.loop, 'subject')
        loop = '';
        par  = 0;
    else
        loop = opt.split.loop;
        par  = opt.split.par;
    end
    
    % --- Detect noise model
    if isfield(dat, 'model')
        noisemodel = dat.model;
    else
        noisemodel = opt.model;
    end
    
    % -------------------------
    % CASE 1: observed velocity
    % -------------------------
    if dat.v.observed
        
        % Lower bound
        % -----------
        % Log-likelihood of a multivariate normal distribution
        % mean = W*z, precision = l*L
        K = prod(opt.tpl.lat)*3;
        if opt.pg.provided
            wz = reconstructVelocity('latent', dat.z.z, 'subspace', model.pg.w, ...
                                     'loop', loop, 'par', par);
            r = numeric(dat.v.v) - wz;
            clear wz
            m = spm_diffeo('vel2mom', r, double([opt.tpl.vs opt.pg.prm]));
            dat.v.lb.reg = r(:)' * m(:);
            clear r m
        else
            dat.v.lb.reg   = single(dat.v.v(:))' * single(dat.v.m(:));
        end
        dat.v.lb.uncty = trace(dat.z.S * model.pg.ww);
        dat.v.lb.val   = -0.5*( K*log(2*pi/model.v.l) - opt.pg.ld ...
                                + model.v.l * dat.v.lb.reg ...
                                + model.v.l * dat.v.lb.uncty );
        dat.v.lb.type  = 'll';
    
    % -----------------------
    % CASE 2: latent velocity
    % -----------------------
    else
        
        % Creator function
        % ----------------
        switch lower(mode)
            case 'zero'
                create = @zeros;
            case 'rand'
                create = @randn;
            otherwise
                error('Unknwon mode %s', mode)
        end
        
        % Allocate and initialise
        % -----------------------
        dat.v.ok = 1; % for GN failure tracking
        dat.v.ok2 = 0; % for GN failure tracking
        dat.v.v = prepareOnDisk(dat.v.v, [opt.tpl.lat 3]);
        dat.v.m = prepareOnDisk(dat.v.m, [opt.tpl.lat 3]);
        if strcmpi(mode, 'zero')
            dat.v.v(:)   = 0;
            dat.v.m(:)   = 0;
            dat.v.lb.reg = 0;
        else
            v = create([opt.tpl.lat 3], 'single');
            dat.v.v(:,:,:,:) = v;
            m = spm_diffeo('vel2mom', v, double([opt.tpl.vs opt.pg.prm]));
            dat.v.m(:,:,:,:) = m;
            if opt.pg.provided
                clear m
                wz = reconstructVelocity('latent', dat.z.z, 'subspace', model.pg.w, ...
                                         'loop', loop, 'par', par);
                r = v - wz;
                clear v wz
                m = spm_diffeo('vel2mom', r, double([opt.tpl.vs opt.pg.prm]));
                dat.v.lb.reg = r(:)' * m(:);
                clear r m
            else
                dat.v.lb.reg = v(:)' * m(:);
                clear v m
            end
        end
        dat.v.lb.uncty = trace(dat.z.S * model.pg.ww);
        
        % Hessian
        % -------
        bx = dat.f.bb.x;
        by = dat.f.bb.y;
        bz = dat.f.bb.z;
        h = zeros([opt.tpl.lat 6], 'single');
        h(bx,by,bz,:) = ghMatchingVel(noisemodel, ...
            model.tpl.mu, dat.f.pf, dat.f.c, model.tpl.gmu, ...
            'bb', dat.f.bb, 'hessian', true, ...
            'loop', loop, 'par', par, 'debug', opt.ui.debug);
        
        % Lower bound
        % -----------
        % KL divergence between multivariate normal distributions
        % posterior: mean = v,    precision = H + l*L
        % prior:     mean = W*z , precision = l*L
        K = prod(opt.tpl.lat)*3;
        dat.v.lb.tr = spm_diffeo('trapprox', h, double([opt.tpl.vs model.v.l*opt.pg.prm]));
        dat.v.lb.tr = dat.v.lb.tr(1);
        dat.v.lb.tr = dat.v.lb.tr / model.v.l;
        clear h
        dat.v.lb.val = -0.5*( model.v.l * dat.v.lb.tr ...
                              - K*log(model.v.l) ...
                              - opt.pg.ld ...
                              + model.v.l * dat.v.lb.reg ...
                              + model.v.l * dat.v.lb.uncty ...
                              - K );
        dat.v.lb.type = 'kl';
        % NOTE: it's missing the term -logdet(H + l*L), but i do not know
        % how to compute it.
        
    end

end

function [dat, model] = batchInitVelocity(mode, dat, model, opt)

    % --- Detect parallelisation scheme
    if strcmpi(opt.split.loop, 'subject') && opt.split.par > 0
        batch = opt.split.batch;
    else
        batch = 1;
    end

    % Init log-likelihood
    % -------------------
    model.lb.v1.val = 0;
    model.lb.v2.val = 0;
    
    % --- Batch processing
    N = numel(dat);
    if opt.ui.verbose, before = plotBatchBegin('Init Vel'); end
    for i=1:ceil(N/batch)
        n1 = (i-1)*batch + 1;
        ne = min(N, i*batch);

        if opt.ui.verbose, before = plotBatch(i, batch, N, 50, before); end
    
        dat(n1:ne) = distribute('OneInitVelocity', dat(n1:ne), model, opt, mode);
        
        for n=n1:ne
            
            % Add individual contributions
            % ----------------------------
            if dat(n).v.observed
                model.lb.v2.val = model.lb.v2.val + dat(n).v.lb.val;
            else
                model.lb.v1.val = model.lb.v1.val + dat(n).v.lb.val;
            end
        end
        
    end
    if opt.ui.verbose, plotBatchEnd; end
end

% -------------------------------------------------------------------------
%    Init Latent
% -------------------------------------------------------------------------

function [dat, model] = batchInitLatent(mode, dat, model, opt)

    % --- Don't parallelise here, there is no need
    
    if strcmpi(mode, 'auto')
        if opt.pg.provided
            mode = 'zero';
        else
            mode = 'rand';
        end
    end
    switch lower(mode)
        case 'zero'
            create = @zeros;
        case 'rand'
            create = @randn;
        otherwise
            error('Unknown mode %s', mode)
    end
    
    % Init model
    % ----------
    N              = numel(dat);
    model.z.z      = zeros(opt.pg.K, 1);
    model.z.zz     = zeros(opt.pg.K);
    model.z.S      = zeros(opt.pg.K);
    model.z.Z      = zeros(opt.pg.K, opt.f.N + opt.v.N);
    
    % Init subjects
    % -------------
    if opt.ui.verbose, before = plotBatchBegin('Init Z'); end
    for n=1:N
        if opt.ui.verbose, before = plotBatch(n, 1, N, 50, before); end
        z = create([opt.pg.K, 1]);
        dat(n).z.z  = z;
        dat(n).z.zz = z*z';
        dat(n).z.S = inv(model.v.l * model.pg.ww + model.z.A);
        model.z.z  = model.z.z  + dat(n).z.z;
        model.z.zz = model.z.zz + dat(n).z.zz;
        model.z.S  = model.z.S  + dat(n).z.S;
        model.z.Z(:,n) = dat(n).z.z;
    end
    if opt.ui.verbose, plotBatchEnd; end
    
    if ~strcmpi(mode, 'zero')
        
        % Center subjects
        % ---------------
        if opt.ui.verbose, before = plotBatchBegin('Center Z'); end
        for n=1:N
            if opt.ui.verbose, before = plotBatch(n, 1, N, 50, before); end
            z = dat(n).z.z - model.z.z/N;
            dat(n).z.z  = z;
            dat(n).z.zz = z*z';
            model.z.Z(:,n) = dat(n).z.z;
        end
        if opt.ui.verbose, plotBatchEnd; end
        model.z.zz = model.z.zz - model.z.z*model.z.z'/N;
        model.z.z  = zeros(opt.pg.K, 1);

        % Orthogonalise
        % -------------
        [U,S] = svd(model.z.zz);
        Rz    = 0.1*sqrt(N/opt.pg.K)*U/diag(sqrt(diag(S)+eps));
        model.z.zz = Rz' * model.z.zz * Rz;
        if opt.pg.provided
            model.z.S  = Rz' * model.z.S  * Rz;
        end
    
        % Center subjects
        % ---------------
        if opt.ui.verbose, before = plotBatchBegin('Rotate Z'); end
        for n=1:N
            if opt.ui.verbose, before = plotBatch(n, 1, N, 50, before); end
            z = Rz' * dat(n).z.z;
            dat(n).z.z  = z;
            dat(n).z.zz = z*z';
            if opt.pg.provided
                dat(n).z.S  = Rz' * dat(n).z.S * Rz;
            else
                dat(n).z.S = inv(model.z.A);
            end
            model.z.Z(:,n) = dat(n).z.z;
        end
        if opt.ui.verbose, plotBatchEnd; end
        
    end
    
    % Lower bound
    % -----------
    model.lb.z.val = 0;
    if opt.ui.verbose, before = plotBatchBegin('LB Z'); end
    for n=1:N
        if opt.ui.verbose, before = plotBatch(n, 1, N, 50, before); end
        dat(n).z.lb.type = 'kl';
        dat(n).z.lb.val  = -0.5*( trace((dat(n).z.S + dat(n).z.zz) * model.z.A) ...
                                  - spm_matcomp('LogDet', model.z.A) ...
                                  - spm_matcomp('LogDet', dat(n).z.S) ...
                                  - opt.pg.K );
        model.lb.z.val = model.lb.z.val + dat(n).z.lb.val;
    end
    if opt.ui.verbose, plotBatchEnd; end

end

% -------------------------------------------------------------------------
%    Rotate
% -------------------------------------------------------------------------

function [dat, model] = batchRotateSubspace(R, iR, dat, model, opt)

    % --- Don't parallelise here, there is no need
    
    % Rotate W
    % --------
    for z=1:size(model.pg.w,3)
        w1  = single(model.pg.w(:,:,z,:,:));
        dim = [size(w1) 1 1 1];
        w1  = reshape(w1, [], opt.pg.K) * iR;
        model.pg.w(:,:,z,:,:) = reshape(w1, dim);
    end
    
    % Rotate sufficient statistics
    model.pg.ww = iR' * model.pg.ww * iR;
    model.z.z   = R   * model.z.z;
    model.z.zz  = R   * model.z.zz  * R';
    model.z.S   = R   * model.z.S   * R';
    model.z.Z   = R   * model.z.Z;
    
    % Rotate subjects
    % ---------------
    N = numel(dat);
    if opt.ui.verbose, before = plotBatchBegin('Rotate Z'); end
    for n=1:N
        if opt.ui.verbose, before = plotBatch(n, 1, N, 50, before); end
        z = R * dat(n).z.z;
        dat(n).z.z  = z;
        dat(n).z.zz = z*z';
        dat(n).z.S  = R * dat(n).z.S * R';
    end
    if opt.ui.verbose, plotBatchEnd; end
end

% -------------------------------------------------------------------------
%    E-step
% -------------------------------------------------------------------------

% ------
% Affine
% ------

function dat = oneFitAffine(dat, model, opt)

    if ~dat.f.observed
        return
    end

    % Detect parallelisation scheme
    % -----------------------------
    if strcmpi(opt.split.loop, 'subject')
        loop = '';
        par  = 0;
        verbose = false;
    else
        loop    = opt.split.lopp;
        par     = opt.split.par;
        verbose = opt.split.verbose;
    end
    
    if isfield(dat, 'model')
        noisemodel = dat.model;
    else
        noisemodel = opt.model;
    end
    
    if opt.q.Mr || dat.q.ok >= 0
        % Compute phi/jac (needed for Affine fitting)
        % ---------------
        if isfield(dat, 'v')
            [dat.v.iphi, dat.v.phi, dat.v.jac] = exponentiateVelocity(dat.v.v, ...
                'iphi', 'phi', 'jac', ...
                'itgr', opt.iter.itg, 'vs', opt.tpl.vs, ...
                'prm', opt.pg.prm,   'bnd', opt.pg.bnd, ...
                'debug', opt.ui.debug, ...
                'output', {dat.v.iphi, dat.v.phi, dat.v.jac});
            iphi = dat.v.iphi;
            phi  = dat.v.phi;
            jac  = dat.v.jac;
        else
            iphi = warps('identity', opt.lat);
            phi  = [];
            jac  = [];
        end
    end
        
    % Penalise previous failure
    % -------------------------
    if dat.q.ok < 0
        dat.q.ok = dat.q.ok + 1;
    else

        % Gauss-Newton iterations
        % -----------------------
        % It is useful to actually find a mode of the posterior (and not only
        % an improved value) when we use the Laplace precision for the update
        % of W. In that case, setting gnit > 1 might help converge faster.
        cumok = false;
        for i=1:opt.iter.gn

            % Compute gradient/hessian
            % ------------------------

            [dat.q.g, dat.q.h] = ghMatchingAffine(noisemodel, ...
                model.tpl.mu, dat.f.pf, dat.f.c, ...
                model.tpl.gmu, dat.q.A, opt.q.B, ...
                phi, jac, 'bb', dat.f.bb, ...
                'Mmu', model.tpl.M, 'loop', loop, 'par', par, ...
                'debug', opt.ui.debug, 'approx', opt.q.hapx);

            if checkarray(model.q.A)
                rind = opt.q.rind;
                [gq, hq] = ghPriorAffine(dat.q.q(rind), model.q.A, 'debug', opt.debug);
                dat.q.g(rind)      = dat.q.g(rind)      + gq;
                dat.q.h(rind,rind) = dat.q.h(rind,rind) + hq;
                clear gq hq
            else
                rind = [];
            end

            dat.q.h = loadDiag(dat.q.h); % Additional regularisation for robustness)

            % Compute search direction
            % ------------------------
            dq = - dat.q.h \ dat.q.g;

            % Line search
            % -----------
            [okq, q, llm, ~, A, pf, c, bb] = lsAffine(...
                noisemodel, dq, dat.q.q, dat.f.lb.val, model.tpl.mu, dat.f.f, ...
                'B', opt.q.B, 'regq', model.q.A, 'rind', rind, ...
                'iphi', iphi, 'Mf', dat.f.M, 'Mmu', model.tpl.M, ...
                'nit', opt.iter.ls, 'loop', loop, 'par', par, ...
                'verbose', verbose, 'debug', opt.ui.debug);

            % Store better values
            % -------------------
            cumok = cumok || okq;
            compute_hessian = okq;
            if okq
                dat.q.q       = q;
                dat.q.qq      = q * q';
                dat.f.lb.val  = llm;
                dat.q.A       = A;
                dat.f.pf      = prepareOnDisk(dat.f.pf, size(pf));
                dat.f.pf(:)   = pf(:);
                dat.f.c       = prepareOnDisk(dat.f.c, size(c));
                dat.f.c(:)    = c(:);
                dat.f.bb      = bb;
            else
                break
            end

        end % < GN iterations
        if cumok
            dat.q.ok2 = 0;
            dat.q.ok  = 1; 
        else
            dat.q.ok2 = dat.q.ok2 - 1;
            dat.q.ok  = dat.q.ok2; 
        end
        
        if opt.q.Mr
            
            % Prior / KL-divergence
            % ---------------------

            if compute_hessian
                dat.q.h = ghMatchingAffine(noisemodel, ...
                    model.tpl.mu, dat.f.pf, dat.f.c, ...
                    model.tpl.gmu, dat.q.A, opt.q.B, ...
                    dat.v.phi, dat.v.jac, ...
                    'bb', dat.f.bb', 'hessian', true, ...
                    'Mmu', model.tpl.M, 'loop', loop, 'par', par, ...
                    'debug', opt.ui.debug, 'approx', opt.q.hapx);

                if checkarray(model.q.A)
                    [~, hq] = ghPriorAffine(dat.q(rind), model.q.A, 'debug', opt.ui.debug);
                    dat.q.h(rind,rind) = dat.q.h(rind,rind) + hq;
                    clear hq
                end

                dat.q.h = loadDiag(dat.q.h); % Additional regularisation for robustness
            end
            dat.q.S = inv(dat.q.h);
        else
            dat.q.S = 0;
        end
    
    end % < penalise previous failure
    
    % Cleaning
    % --------
    % I should probably clear variables and remove files that are not
    % useful anymore. This will cause less disk and broadband usage.
    dat.q.g    = rmarray(dat.q.g);
    dat.q.h    = rmarray(dat.q.h);
    dat.v.iphi = rmarray(dat.v.ipsi);
    dat.v.ipsi = rmarray(dat.v.iphi);
    dat.v.phi  = rmarray(dat.v.phi);
    dat.v.jac  = rmarray(dat.v.jac);
end

function [dat, model] = batchFitAffine(dat, model, opt)

    % --- Detect parallelisation scheme
    if strcmpi(opt.split.loop, 'subject') && opt.split.par > 0
        batch = opt.split.batch;
    else
        batch = 1;
    end

    % Init gradient/hessian
    % ---------------------
    model.lb.m.val = 0;
    model.q.q      = zeros(opt.q.M, 1);
    model.q.qq     = zeros(opt.q.M);
    model.q.S      = zeros(opt.q.M);
    
    N   = numel(dat);
    Nok = 0;
    for n=1:N
        if isfield(dat(n).q, 'ok')
            Nok = Nok + (dat(n).q.ok >= 0);
        end
    end
    
    % --- Batch processing
    if opt.ui.verbose, before = plotBatchBegin('Fit Affine'); end
    for i=1:ceil(N/batch)
        n1 = (i-1)*batch + 1;
        ne = min(N, i*batch);

        if opt.ui.verbose, before = plotBatch(i, batch, N, 50, before); end
    
        % Compute subjects grad/hess w.r.t. initial velocity
        % --------------------------------------------------
        dat(n1:ne) = distribute('oneFitAffine', dat(n1:ne), model, opt);
        
        for n=n1:ne
            
            if dat(n).f.observed
                % Add individual contributions
                % ----------------------------
                model.lb.m.val = model.lb.m.val + dat(n).f.lb.val;
                model.q.q      = model.q.q  + dat(n).q.q;
                model.q.qq     = model.q.qq + dat(n).q.qq;
                model.q.S      = model.q.S  + dat(n).q.S;
            end
        end
        
    end
    okq = 0;
    for n=1:N
        if isfield(dat(n).q, 'ok')
            okq = okq + (dat(n).q.ok > 0);
        end
    end
    if opt.ui.verbose, fprintf(' | %3d / %3d / %3d', okq, Nok, N); end
    if opt.ui.verbose, plotBatchEnd; end

end

% ------------------
% Latent coordinates
% ------------------

function dat = oneFitLatent(dat, model, opt)
    
    % Update covariance
    % -----------------
    dat.z.S = model.v.l * model.pg.ww + model.z.A;
    dat.z.S = spm_matcomp('Inv', dat.z.S);
    
    % Update mean
    % -----------
    wm = zeros([opt.pg.K 1]);
    m  = numeric(dat.v.m);
    for k=1:opt.pg.K
        w1    = model.pg.w(:,:,:,:,k);
        wm(k) = w1(:)'*m(:);
    end
    clear m w1
    dat.z.z  = model.v.l * dat.z.S * wm;
    dat.z.zz = dat.z.z * dat.z.z';
    
    % Lower bound
    % -----------
    dat.z.lb.val  = -0.5*( trace((dat.z.S + dat.z.zz) * model.z.A) ...
                           - spm_matcomp('LogDet', model.z.A) ...
                           - spm_matcomp('LogDet', dat.z.S) ...
                           - opt.pg.K );
    
end

function [dat, model] = batchFitLatent(dat, model, opt)

    % --- Detect parallelisation scheme
    if strcmpi(opt.split.loop, 'subject') && opt.split.par > 0
        batch = opt.split.batch;
    else
        batch = 1;
    end
    N = numel(dat);

    % Init gradient/hessian
    % ---------------------
    model.z.lb.val = 0;
    model.z.z      = zeros(opt.pg.K, 1);
    model.z.zz     = zeros(opt.pg.K);
    model.z.S      = zeros(opt.pg.K);
    model.z.Z      = zeros(opt.pg.K, opt.f.N+opt.v.N);
    
    % --- Batch processing
    if opt.ui.verbose, before = plotBatchBegin('Fit Latent'); end
    for i=1:ceil(N/batch)
        n1 = (i-1)*batch + 1;
        ne = min(N, i*batch);

        if opt.ui.verbose, before = plotBatch(i, batch, N, 50, before); end
    
        % Compute subjects grad/hess w.r.t. initial velocity
        % --------------------------------------------------
        dat(n1:ne) = distribute('oneFitLatent', dat(n1:ne), model, opt);
        
        for n=n1:ne
            
            % Add individual contributions
            % ----------------------------
            model.z.z      = model.z.z      + dat(n).z.z;
            model.z.zz     = model.z.zz     + dat(n).z.zz;
            model.z.S      = model.z.S      + dat(n).z.S;
            model.lb.z.val = model.lb.z.val + dat(n).z.lb.val;
            model.z.Z(:,n) = dat(n).z.z;
            
        end
        
    end
    if opt.ui.verbose, plotBatchEnd; end
end

% --------------
% Velocity field
% --------------

function dat = oneFitVelocity(dat, model, opt)

    % Detect parallelisation scheme
    % -----------------------------
    if strcmpi(opt.split.loop, 'subject')
        loop    = '';
        par     = 0;
        verbose = false;
    else
        loop    = opt.split.loop;
        par     = opt.split.par;
        verbose = opt.ui.verbose;
    end
    
    % If v is observed
    % ----------------
    if dat.v.observed
        return
    end
    
    if isfield(dat, 'model')
        noisemodel = dat.model;
    else
        noisemodel = opt.model;
    end
    
    if isfield(dat, 'q') && isfield(dat.q, 'A'),  A = dat.q.A;
    else,                                         A = eye(4);  end
    
    % Set L boundary conditions
    % -------------------------
    spm_diffeo('boundary', opt.pg.bnd);
        
    % Penalise previous failure
    % -------------------------
    if dat.v.ok < 0
        dat.v.ok = dat.v.ok + 1;
    else
        
        % Gauss-Newton iterations
        % -----------------------
        cumok = false;
        for i=1:opt.iter.gn

            % Compute gradient/hessian
            % ------------------------
            g = zeros([opt.tpl.lat 3], 'single'); 
            h = zeros([opt.tpl.lat 6], 'single'); 
            bx = dat.f.bb.x;
            by = dat.f.bb.y;
            bz = dat.f.bb.z;

            [g(bx,by,bz,:), h(bx,by,bz,:)] = ghMatchingVel(...
                noisemodel, ...
                model.tpl.mu, dat.f.pf, dat.f.c, model.tpl.gmu, ...
                'bb', dat.f.bb, 'loop', loop, 'par', par, 'debug', opt.ui.debug);

            wz = reconstructVelocity('latent', dat.z.z, 'subspace', model.pg.w, ...
                                     'loop', loop, 'par', par);
            v = numeric(dat.v.v);
            r = v - wz;
            clear wz
            g = g + ghPriorVel(r, opt.tpl.vs, model.v.l * opt.pg.prm, opt.pg.bnd);

            % Compute search direction
            % ------------------------
            dv = -spm_diffeo('fmg', single(h), single(g), ...
                double([opt.tpl.vs model.v.l * opt.pg.prm 2 2]));
            clear g

            % Line search
            % -----------
            result = lsVelocity2(...
                noisemodel, dv, r, dat.f.lb.val, model.tpl.mu, dat.f.f, ...
                'v0', v, 'lam', model.v.l, 'prm', opt.pg.prm, ...
                'itgr', opt.iter.itg, 'bnd', opt.pg.bnd, ...
                'A', A, 'Mf', dat.f.M, 'Mmu', model.tpl.M, ...
                'match', 'push', ...
                'nit', opt.iter.ls,  'par', par, 'loop', loop, ...
                'verbose', verbose, 'debug', opt.ui.debug, ...
                'pf', dat.f.pf, 'c', dat.f.c, 'wa', dat.tpl.wa, 'wmu', dat.tpl.wmu);
%                 'match', 'pull', 'itrp', opt.tpl.itrp, 'tplbnd', opt.tpl.bnd, ...

            % Store better values
            % -------------------
            cumok = cumok || result.ok;
            compute_hessian = result.ok;
            if result.ok
                dat.f.lb.val  = result.match;
                dat.v.v       = prepareOnDisk(dat.v.v, size(result.v));
                dat.v.v(:)    = result.v(:);
                m = spm_diffeo('vel2mom', result.v, double([opt.tpl.vs opt.pg.prm]));
                dat.v.m       = prepareOnDisk(dat.v.m, size(m));
                dat.v.m(:)    = m(:);
                clear m
                if isfield(result, 'pf')
                    dat.f.pf = result.pf;
                end
                if isfield(result, 'c')
                    dat.f.c = result.c;
                end
                if isfield(result, 'bb')
                    dat.f.bb = result.bb;
                end
                if isfield(result, 'wa')
                    dat.tpl.wmu = result.wmu;
                end
                if isfield(result, 'wmu')
                    dat.tpl.wmu = result.wmu;
                end
                r = result.r;
            else
                break
            end
            clear result

        end % < GN iterations
        if cumok
            dat.v.ok2 = 0;
            dat.v.ok  = 1; 
        else
            dat.v.ok2 = dat.v.ok2 - 1;
            dat.v.ok  = dat.v.ok2; 
        end
    

        % -----------
        % Lower bound
        % -----------
        % KL divergence between multivariate normal distributions
        % posterior: mean = v,    precision = H + l*L
        % prior:     mean = W*z , precision = l*L
        
        % Regularisation part
        % -------------------
        m = spm_diffeo('vel2mom', r, double([opt.tpl.vs opt.pg.prm]));
        dat.v.lb.reg = r(:)' * m(:);
        clear r m
        
        % Trace(P\L) part
        % ---------------
        if compute_hessian
            bx = dat.f.bb.x;
            by = dat.f.bb.y;
            bz = dat.f.bb.z;
            h  = zeros([opt.tpl.lat 6], 'single'); 
            h(bx,by,bz,:) = ghMatchingVel(noisemodel, ...
                model.tpl.mu, dat.f.pf, dat.f.c, model.tpl.gmu, ...
                'bb', dat.f.bb, 'hessian', true, ...
                'loop', loop, 'par', par, 'debug', opt.ui.debug);
        end
        dat.v.lb.tr = spm_diffeo('trapprox', h, double([opt.tpl.vs model.v.l * opt.pg.prm]));
        dat.v.lb.tr = dat.v.lb.tr(1);
        dat.v.lb.tr = dat.v.lb.tr / model.v.l;
        
        % LogDet(P)
        % ---------
        % Approximation where all off-diagonal elements of L are zero
        h(:,:,:,1) = h(:,:,:,1) * model.v.l * opt.pg.ker(1);
        h(:,:,:,2) = h(:,:,:,2) * model.v.l * opt.pg.ker(2);
        h(:,:,:,3) = h(:,:,:,3) * model.v.l * opt.pg.ker(3);
        h = spm_matcomp('Pointwise3', h, 'd');
        h(h <= 0) = nan;
        dat.v.lb.ld = sum(log(h(:)), 'omitnan');
        clear h
        
        % KL divergence
        % -------------
        K = prod(opt.tpl.lat) * 3;
        dat.v.lb.val = -0.5*( model.v.l * dat.v.lb.tr ...
                              - K*log(model.v.l) ...
                              - opt.pg.ld ...
                              + dat.v.lb.ld ...
                              + model.v.l * dat.v.lb.reg ...
                              + model.v.l * dat.v.lb.uncty ...
                              - K );
    
    end % < penalise previous failure
    
    % Cleaning
    % --------
    % Just in case
    dat.v.iphi = rmarray(dat.v.iphi);
    dat.v.ipsi = rmarray(dat.v.ipsi);
    dat.tpl.wa = rmarray(dat.tpl.wa);
end

function [dat, model] = batchFitVelocity(dat, model, opt)

    % --- Detect parallelisation scheme
    if strcmpi(opt.split.loop, 'subject') && opt.split.par > 0
        batch = opt.split.batch;
    else
        batch = 1;
    end

    % Init gradient/hessian
    % ---------------------
    model.lb.m.val  = 0;
    model.lb.v1.val = 0;
    model.lb.v2.val = 0;
    model.v.tr      = 0;
    model.v.reg     = 0;
    model.v.uncty   = 0;
    
    total = opt.f.N;
    N     = numel(dat);
    okpre = 0;
    for n=1:N
        if isfield(dat(n).v, 'ok')
            okpre = okpre + (dat(n).v.ok >= 0);
        end
    end
    
    % --- Batch processing
    if opt.ui.verbose, before = plotBatchBegin('Fit Vel'); end
    for i=1:ceil(N/batch)
        n1 = (i-1)*batch + 1;
        ne = min(N, i*batch);

        if opt.ui.verbose, before = plotBatch(i, batch, N, 50, before); end
    
        % Compute subjects grad/hess w.r.t. initial velocity
        % --------------------------------------------------
        dat(n1:ne) = distribute('oneFitVelocity', dat(n1:ne), model, opt);
        
        
        for n=n1:ne
            
            % Add individual contributions
            % ----------------------------
            model.v.reg     = model.v.reg    + dat(n).v.lb.reg;
            model.v.uncty   = model.v.uncty  + dat(n).v.lb.uncty;
            if dat(n).v.observed
                model.lb.v2.val = model.lb.v2.val + dat(n).v.lb.val;
            else
                model.lb.m.val  = model.lb.m.val  + dat(n).f.lb.val;
                model.lb.v1.val = model.lb.v1.val + dat(n).v.lb.val;
                model.v.tr      = model.v.tr      + dat(n).v.lb.tr;
            end
            
        end
        
    end
    okpost = 0;
    for n=1:N
        if isfield(dat(n).v, 'ok')
            okpost = okpost + (dat(n).v.ok > 0);
        end
    end
    if opt.ui.verbose, fprintf(' | %3d / %3d / %3d', okpost, okpre, total); end
    if opt.ui.verbose, plotBatchEnd; end

end

% -------------------------------------------------------------------------
%    Update Lower Bound
% -------------------------------------------------------------------------

function dat = oneLB(dat, model, opt, which)
    
    % Detect parallelisation scheme
    % -----------------------------
    if strcmpi(opt.split.loop, 'subject')
        loop    = '';
        par     = 0;
    else
        loop    = opt.split.loop;
        par     = opt.split.par;
    end
    if isfield(dat, 'model')
        noisemodel = dat.model;
    else
        noisemodel = opt.model;
    end
        
    
    if strcmpi(which, 'matching')
        
        % Update matching terms
        % ---------------------
        if dat.f.observed
%             dat.f.lb.val = llMatching(noisemodel, dat.tpl.wmu, dat.f.f, ...
%                                       'par', par, 'loop', loop, ...
%                                       'debug', opt.ui.debug);
            dat.f.lb.val = llMatching(noisemodel, model.tpl.mu, dat.f.pf, dat.f.c, ...
                                      'par', par, 'loop', loop, ...
                                      'debug', opt.ui.debug, 'bb', dat.f.bb);
        end
        return
        
    end
    
    if strcmpi(which, 'precisionz')
        
        dat.z.lb.val = -0.5*( trace((dat.z.zz + dat.z.S)*model.z.A) ...
                              - spm_matcomp('LogDet', model.z.A) ...
                              - spm_matcomp('LogDet', dat.z.S) ...
                              - opt.pg.K );
        return
        
    end
    
    if strcmpi(which, 'precisionq')
        
        if ~opt.f.observed
            return
        end
        rind = opt.q.rind;
        dat.q.lb.val = -0.5*( trace((dat.q.qq(rind,rind) + dat.q.S(rind,rind))*model.q.A) ...
                              - spm_matcomp('LogDet', model.q.A) ...
                              - spm_matcomp('LogDet', dat.q.S(rind,rind)) ...
                              - opt.q.Mr );
        return
        
    end
    
    if strcmpi(which, 'orthogonalise')
        
        % Update after orthogonalisation
        % ------------------------------
        % Here, no need to recompute terms depending on Wz
        dat.v.lb.uncty = trace(dat.z.S * model.pg.ww);
        dat.z.lb.val  = -0.5*( trace((dat.z.S + dat.z.zz) * model.z.A) ...
                               - spm_matcomp('LogDet', model.z.A) ...
                               - spm_matcomp('LogDet', dat.z.S) ...
                               - opt.pg.K );
        

    elseif any(strcmpi(which, {'latent', 'subspace'}))
        
        % Update regularisation and uncertainty terms
        % -------------------------------------------
        wz = reconstructVelocity('latent', dat.z.z, 'subspace', model.pg.w, ...
                                 'loop', loop, 'par', par);
        r = numeric(dat.v.v) - wz;
        clear wz
        m = spm_diffeo('vel2mom', r, double([opt.tpl.vs opt.pg.prm]));
        dat.v.lb.reg = r(:)' * m(:);
        clear r m
        dat.v.lb.uncty = trace(dat.z.S * model.pg.ww);
    end
    
    % Update complete LL or KL
    % ------------------------
    % cases latent/subspace/orthogonalise/lambda
    K = prod(opt.tpl.lat) * 3;
    if dat.v.observed
        dat.v.lb.val = -0.5*( K*log(2*pi/model.v.l) ...
                              - opt.pg.ld ...
                              + model.v.l * dat.v.lb.reg ...
                              + model.v.l * dat.v.lb.uncty );
    else
        dat.v.lb.val = -0.5*( - K - K*log(model.v.l) ...
                              - opt.pg.ld + dat.v.lb.ld ...
                              + model.v.l * dat.v.lb.reg ...
                              + model.v.l * dat.v.lb.uncty ...
                              + model.v.l * dat.v.lb.tr );
    end

end

function [dat, model] = batchLB(which, dat, model, opt)

    % --- Detect parallelisation scheme
    if strcmpi(opt.split.loop, 'subject') && opt.split.par > 0
        batch = opt.split.batch;
    else
        batch = 1;
    end

    % Init gradient/hessian
    % ---------------------
    model.lb.m.val  = 0;
    model.lb.z.val  = 0;
    model.lb.q.val  = 0;
    model.lb.v1.val = 0;
    model.lb.v2.val = 0;
    model.v.tr      = 0;
    model.v.reg     = 0;
    model.v.uncty   = 0;
    
    N = numel(dat);
    
    % --- Batch processing
    if opt.ui.verbose, before = plotBatchBegin(['LB ' which(1:min(numel(which),3))]); end
    for i=1:ceil(N/batch)
        n1 = (i-1)*batch + 1;
        ne = min(N, i*batch);

        if opt.ui.verbose, before = plotBatch(i, batch, N, 50, before); end
    
        % Compute subjects grad/hess w.r.t. initial velocity
        % --------------------------------------------------
        dat(n1:ne) = distribute('oneLB', dat(n1:ne), model, opt, which);
        
        for n=n1:ne
            
            % Add individual contributions
            % ----------------------------
            model.lb.z.val  = model.lb.z.val + dat(n).z.lb.val;
            model.v.reg     = model.v.reg    + dat(n).v.lb.reg;
            model.v.uncty   = model.v.uncty  + dat(n).v.lb.uncty;
            if dat(n).v.observed
                model.lb.v2.val = model.lb.v2.val + dat(n).v.lb.val;
            else % dat(n).f.observed
                model.lb.m.val  = model.lb.m.val  + dat(n).f.lb.val;
                model.lb.v1.val = model.lb.v1.val + dat(n).v.lb.val;
                model.v.tr      = model.v.tr      + dat(n).v.lb.tr;
                model.lb.q.val  = model.lb.q.val + dat(n).q.lb.val;
            end
            
        end
        
    end
    if opt.ui.verbose, plotBatchEnd; end

    % Model specific parts
    % --------------------
    switch lower(which)
        case {'subspace', 'orthogonalise'}
            model.lb.w.val = llPriorSubspace(model.pg.w, model.pg.ww, opt.pg.ld);
        case 'precisionz'
            if opt.z.n0
                model.lb.Az.val = -spm_prob('Wishart', 'kl', ...
                                           model.z.A,   opt.z.n0+opt.v.N+opt.f.N, ...
                                           opt.z.A0,    opt.z.n0, ...
                                           'normal');
            end
        case 'precisionq'
            if opt.q.n0 && opt.q.Mr
                model.lb.Aq.val = -spm_prob('Wishart', 'kl', ...
                                           model.q.A,   opt.q.n0+opt.f.N, ...
                                           opt.q.A0,    opt.q.n0, ...
                                           'normal');
            end
        case 'lambda'
            if opt.v.n0
                model.lb.l.val  = -spm_prob('Gamma', 'kl', ...
                                            model.v.l, opt.v.N+opt.f.N+opt.v.n0, ...
                                            opt.v.l0,  opt.v.n0, ...
                                            prod(opt.tpl.lat)*3, 'normal');
            end
    end
end


% =========================================================================

function dat = distribute(func, dat, model, opt, varargin)

    % --- Convert function name to function handle
    if ischar(func)
        funcname = func;
        func = @(varargin) pgva_batch(funcname, varargin{:});
        funcstr = ['@(varargin) pgva_batch(''' funcname ''', varargin{:})']; 
    else
        funcstr = func2str(func);
    end

    % --- No distribution
    if ~strcmpi(opt.split.loop, 'subject')
        if ~isempty(varargin)
            for n=1:numel(dat)
                dat(n) = func(dat(n), model, opt, varargin{:});
            end
        else
            for n=1:numel(dat)
                dat(n) = func(dat(n), model, opt);
            end
        end
        
    % --- Distribute on cluster
    elseif isfield(opt, 'cluster') && opt.cluster
        % - Prepare
        opt.split.loop    = '';
        opt.split.par     = false;
        opt.ui.verbose = true;
        if ~isfield(opt, 'translate_path')
            translate_path = opt.translate_path;
        else
            translate_path = {};
        end
        if translate_path
            [dat, model] = translatePath(dat, model, translate_path{1}, translate_path{2});
            output_dir = strrep(opt.directory, translate_path{1}, translate_path{2});
        else
            output_dir = opt.directory;
        end
        client_tmpdir = fullfile(opt.directory, 'tmp');
        server_tmpdir = fullfile(output_dir, 'tmp');
        if ~exist(client_tmpdir, 'dir')
            mkdir(client_tmpdir);
        end
        client_model = fullfile(client_tmpdir, 'model.mat');
        server_model = fullfile(server_tmpdir, 'model.mat');
        save(client_model, 'model', 'opt', 'varargin');
        
        % - Distribute jobs
        for n=1:numel(dat)
            client_subjtmpdir = fullfile(client_tmpdir, str(n));
            server_subjtmpdir = fullfile(server_tmpdir, str(n));
            if ~exist(client_subjtmpdir, 'dir')
                mkdir(client_subjtmpdir);
            end
            client_dat = fullfile(client_subjtmpdir, 'dat.mat');
            server_dat = fullfile(server_subjtmpdir, 'dat.mat');
            foo.dat = dat(n);
            save(client_dat, '-struct', 'foo', 'dat');
            job = [...
                '#$ /bin/sh \n' ...
                '#$ -N job' str(n) '\n'...
                '/share/apps/matlab '...
                '-nojvm -nodesktop -nosplash -singleCompThread ' ...
                '-r ' ...
                'addpath(''/data/' opt.cluster.username '/spm''); ' ...
                'addpath(''/data/' opt.cluster.username '/spm/toolbox/Longitudinal''); '
                'addpath(''/data/' opt.cluster.username '/spm/toolbox/Shoot''); ' ...
                'addpath(''/data/' opt.cluster.username '/auxiliary-functions''); ' ...
                'addpath(genpath(''/data/' opt.cluster.username '/shape-toolbox'')); ' ...
                'load(''' server_model '''); ' ...
                'load(''' server_dat '''); ' ...
                'func = ' funcstr '; ' ...
                'dat = func(dat, model, opt, varargin{:}); ' ...
                'save(''' server_dat ''', dat); ' ...
                ];
            client_job = fullfile(client_subjtmpdir, 'job.sh');
            server_job = fullfile(server_subjtmpdir, 'job.sh');
            dlmwrite(client_job, job, 'dlm', '');
            system(['ssh ' opt.cluster.username '@' opt.cluster.ip ...
                    ' ''qsub -l vf=2G -l h_vmem=2G ' server_job '''']);
        end
        
        % - Wait job (will only be submitted when all other are finished)
        client_job = fullfile(client_tmpdir, 'waitjob.sh');
        server_job = fullfile(server_tmpdir, 'waitjob.sh');
        dlmwrite(client_job, '', 'dlm', '');
        endjob = ['ssh ' opt.cluster.username '@' opt.cluster.ip ...
                  '''qsub -l vf=10K -l h_vmem=10K '];
%         endjob = [endjob '-hold_jid '];
        endjob = [endjob '-sync -W depend=afterany'];
        for n=1:numel(dat)
            endjob = [endjob ':job' str(n)];
        end
%         endjob = endjob(1:end-1);
        endjob = [endjob server_job ''''];
        system(endjob);
        
        % - Read output data
        failed = 0;
        for n=1:numel(dat)
            [~, res] = system(['ssh ' opt.cluster.username '@' opt.cluster.ip ...
                                ' qacct -j job' str(n)]);
            exit_status = str2double(regexp(res, '^\s*exit_status\s*:\s*(\d+)'));
            if exit_status > 0
                failed = failed + 1;
            end
            client_dat = fullfile(client_subjtmpdir, 'dat.mat');
            result = load(client_dat, 'dat');
            dat(n) = result.dat;
        end
        if failed > 0
            error('Distribution failure: %d/%d jobs terminated unexplectedly.', failed, numel(dat));
        end
        dat = translatePath(dat, struct, translate_path{2}, translate_path{1});
      
    % --- Distribute locally  
    else
        if ~isempty(varargin)
            if opt.split.par
                parfor (n=1:numel(dat), double(opt.split.par))
                    dat(n) = func(dat(n), model, opt, varargin{:});
                end
            else
                for n=1:numel(dat)
                    dat(n) = func(dat(n), model, opt, varargin{:});
                end
            end
        else
            if opt.split.par
                parfor (n=1:numel(dat), double(opt.split.par))
                    dat(n) = func(dat(n), model, opt);
                end
            else
                for n=1:numel(dat)
                    dat(n) = func(dat(n), model, opt);
                end
            end
        end
    end
        
end

% =========================================================================

function [dat, model] = translatePath(dat, model, onclient, onserver)

    datfields = fieldnames(dat);
    for n=1:numel(N)
       for i=1:numel(datfields)
           f = datfields{i};
           if isa(dat(n).(f), 'file_array')
               dat(n).(f).fname = strrep(dat(n).(f).fname, onclient, onserver);
           end
       end
    end
    
   if ~isempty(model)
       modelfields = fieldnames(model);
       for i=1:numel(modelfields)
           f = modelfields{i};
           if isa(model.(f), 'file_array')
               model.(f).fname = strrep(model.(f).fname, onclient, onserver);
           end
       end
   end
   
end
        