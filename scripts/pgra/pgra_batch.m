function varargout = pgra_batch(id, varargin)
% _________________________________________________________________________
%
%             Batch functions for the PGRA model (2017)
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
% FORMAT dat = pgra_batch('InitPush', dat, model, opt)
%   Initialise warped images (Change of lattice, but no deformation)
%
% FORMAT [dat, model] = pgra_batch('InitPull', dat, model, opt)
%   Initialise warped template (Change of lattice, but no deformation)
%   + matching log-likelihood
%
% FORMAT [dat, model] = pgra_batch('InitAffine', mode, dat, model, opt)
%   Initialise affine coordinates ('zero'/'rand') + lower bound
%
% FORMAT [dat, model] = pgra_batch('InitResidual', mode, dat, model, opt)
%   Initialise reisudal fields ('zero'/'rand') + lower bound
%
% FORMAT [dat, model] = pgra_batch('InitLatent', mode, dat, model, opt)
%   Initialise latent coordinates ('zero'/'rand') + lower bound
%   If random, insure they are zero-centered and orthogonal.
%
% FORMAT [dat, model] = pgra_batch('InitLowerBound', dat, model, opt)
%   Initialise covariance for Laplace approximations + lower bound
%
% ---
% Fit
% ---
%
% FORMAT [dat, model] = pgra_batch('FitAffine', dat, model, opt)
%   Gauss-Newton update of the affine coordinates
%
% FORMAT [dat, model] = pgra_batch('FitLatent', dat, model, opt)
%   Gauss-Newton update of the latent coordinates
%
% FORMAT [dat, model] = pgra_batch('FitResidual', dat, model, opt)
%   Gauss-Newton update of the residual fields
%
% ----------
% Correction
% ----------
% 
% FORMAT dat = pgra_batch('RotateSubspace', R, iR, dat, model, opt)
%   Apply a rotation matrix to all latent coordinates
%
% % FORMAT dat = pgra_batch('CentreLatent', dat, opt, (mean))
% %   Subtract the mean to all latent coordinates
%
% ------
% Update
% ------
%
% FORMAT [dat, model] = pgra_batch('batchGradHessSubspace', dat, model, opt)
%   Compute gradient and Hessian for subspace update
%
% FORMAT [dat, model] = pgra_batch('LB', var, dat, model, opt)
%   Update or compute the part of the lower bound that depends on var.
%   var can be 'Matching'/'Lambda'/'Subspace'.
% _________________________________________________________________________
%
% The following subfunctions should not be used usually. They are only
% needed when distributing jobs on a cluster.
% 
% FORMAT dat = pgra_batch('OneInitPush',     dat, model, opt)
% FORMAT dat = pgra_batch('OneInitPull',     dat, model, opt)
% FORMAT dat = pgra_batch('OneInitResidual', dat, model, opt, mode)
% FORMAT dat = pgra_batch('OneInitLowerBound', dat, model, opt)
% FORMAT dat = pgra_batch('OneFitAffine',    dat, model, opt)
% FORMAT dat = pgra_batch('OneFitLatent',    dat, model, opt)
% FORMAT dat = pgra_batch('OneFitResidual',  dat, model, opt)
% FORMAT dat = pgra_batch('OneGradHessVelocity',  dat, model, opt)
% FORMAT dat = pgra_batch('OneLB',           dat, model, opt, var)
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
        case 'initresidual'
            [varargout{1:nargout}] = batchInitResidual(varargin{:});
        case 'initlowerbound'
            [varargout{1:nargout}] = batchInitLowerBound(varargin{:});
        case 'rotatesubspace'
            [varargout{1:nargout}] = batchRotateSubspace(varargin{:});
        case 'fitaffine'
            [varargout{1:nargout}] = batchFitAffine(varargin{:});
        case 'fitlatent'
            [varargout{1:nargout}] = batchFitLatent(varargin{:});
        case 'fitresidual'
            [varargout{1:nargout}] = batchFitResidual(varargin{:});
        case 'gradhesssubspace'
            [varargout{1:nargout}] = batchGradHessSubspace(varargin{:});
        case 'lb'
            [varargout{1:nargout}] = batchLB(varargin{:});
        % ------
        %  STEP
        % ------
        case 'oneinitpush'
            [varargout{1:nargout}] = oneInitPush(varargin{:});
        case 'oneinitpull'
            [varargout{1:nargout}] = oneInitPull(varargin{:});
        case 'oneinitresidual'
            [varargout{1:nargout}] = oneInitResidual(varargin{:});
        case 'oneinitlowerbound'
            [varargout{1:nargout}] = oneInitLowerBound(varargin{:});
        case 'onefitaffine'
            [varargout{1:nargout}] = oneFitAffine(varargin{:});
        case 'onefitlatent'
            [varargout{1:nargout}] = oneFitLatent(varargin{:});
        case 'onefitresidual'
            [varargout{1:nargout}] = oneFitResidual(varargin{:});
        case 'onegradhessvelocity'
            [varargout{1:nargout}] = oneGradHessVelocity(varargin{:});
        case 'onelb'
            [varargout{1:nargout}] = oneLB(varargin{:});
    end

end

%% ------------------------------------------------------------------------
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

%% ------------------------------------------------------------------------
%    Init Push
% -------------------------------------------------------------------------

function dat = oneInitPush(dat, model, opt)

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
    if opt.pg.provided
        iphi = exponentiateVelocity(dat.v.v, 'iphi', ...
            'itgr', opt.iter.itg, 'vs', opt.tpl.vs, 'prm', opt.pg.prm);
    else
        iphi = spm_warps('identity', opt.tpl.lat);
    end
    if isfield(dat, 'q') && isfield(dat.q, 'A'),  A = dat.q.A;
    else,                                         A = eye(4);  end
    ipsi = reconstructIPsi(A, iphi, 'lat', opt.tpl.lat, ...
                           'Mf',  dat.f.M, 'Mmu', model.tpl.M, ...
                           'debug', opt.ui.debug);
    clear iphi
    [dat.f.pf, dat.f.c, dat.f.bb] = pushImage(ipsi, dat.f.f, opt.tpl.lat, ...
                                              'output', {dat.f.pf, dat.f.c}, ...
                                              'loop', loop, 'par', par, ...
                                              'debug', opt.ui.debug);
    if isa(dat.v.ipsi, 'file_array')
        dat.v.ipsi = copyarray(ipsi, dat.v.ipsi);
    end
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
    
        [opt.dist, dat(n1:ne)] = distribute(opt.dist, 'pgra_batch', 'OneInitPush', 'inplace', dat(n1:ne), model, opt);
        
    end
    if opt.ui.verbose, plotBatchEnd; end
end

%% ------------------------------------------------------------------------
%    Init Pull
% -------------------------------------------------------------------------

function dat = oneInitPull(dat, model, opt)

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
    if opt.tpl.cat
        dat.tpl.wa = pullTemplate(dat.v.ipsi, model.tpl.a, ...
            'par', par, 'output', dat.tpl.wa, 'debug', opt.ui.debug);
        dat.tpl.wmu = reconstructProbaTemplate(dat.tpl.wa, ...
            'loop', loop, 'par', par, 'output', dat.tpl.wmu, ...
            'debug', opt.ui.debug);
        dat.tpl.wa = rmarray(dat.tpl.wa);
    else
        dat.tpl.wmu = pullTemplate(dat.v.ipsi, model.tpl.mu, ...
            'par', par, 'output', dat.tpl.wmu, 'debug', opt.ui.debug);
    end
    
    % Matching likelihood
    % -------------------
    dat.f.lb.type = 'll';
    dat.f.lb.val  = llMatching(noisemodel, dat.tpl.wmu, dat.f.f, ...
        'loop', loop, 'par', par, 'debug', opt.ui.debug);
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
    
        [opt.dist, dat(n1:ne)] = distribute(opt.dist, 'pgra_batch', 'OneInitPull', 'inplace', dat(n1:ne), model, opt);
        
        for n=n1:ne
            model.lb.m.val = model.lb.m.val + dat(n).f.lb.val;
        end
    end
    if opt.ui.verbose, plotBatchEnd; end
end
%% ------------------------------------------------------------------------
%    Init Laplace
% -------------------------------------------------------------------------

function dat = oneInitLowerBound(dat, model, opt)

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
    
    % --------------
    %    Residual
    % --------------
    if opt.optimise.r.r
    
        % Wz uncertainty
        % --------------
        dat.v.lb.uncty = trace(dat.z.S * model.pg.ww);
        
        % Hessian
        % -------
        if strcmpi(opt.match, 'pull')
            h = ghMatchingVel(noisemodel, ...
                dat.tpl.wmu, dat.f.f, model.tpl.gmu, ...
                'ipsi', dat.v.ipsi, 'hessian', true, ...
                'loop', loop, 'par', par, 'debug', opt.ui.debug);
        else
            h = zeros([opt.tpl.lat 6], 'single');
            h(dat.f.bb.x,dat.f.bb.y,dat.f.bb.z,:) = ghMatchingVel(...
                noisemodel, ...
                model.tpl.mu, dat.f.pf, model.tpl.gmu, ...
                'count', dat.f.c, 'bb', dat.f.bb, 'hessian', true, ...
                'loop', loop, 'par', par, 'debug', opt.ui.debug);
        end

        % Trace
        % -----
        dat.v.lb.tr = spm_diffeo('trapprox', h, double([opt.tpl.vs (model.mixreg.w(1)*model.r.l + model.mixreg.w(2))*opt.pg.prm]));
        dat.v.lb.tr = dat.v.lb.tr(1);
        dat.v.lb.tr = dat.v.lb.tr / (model.mixreg.w(1)*model.r.l + model.mixreg.w(2));

        % LogDet(P)
        % ---------
        % Approximation where all off-diagonal elements of L are zero
        h(:,:,:,1) = h(:,:,:,1) * (model.mixreg.w(1)*model.r.l + model.mixreg.w(2)) * opt.pg.ker(1);
        h(:,:,:,2) = h(:,:,:,2) * (model.mixreg.w(1)*model.r.l + model.mixreg.w(2)) * opt.pg.ker(2);
        h(:,:,:,3) = h(:,:,:,3) * (model.mixreg.w(1)*model.r.l + model.mixreg.w(2)) * opt.pg.ker(3);
        if opt.model.dim == 2
            h(:,:,:,3) = 1;
        end
        h = spm_matcomp('Pointwise3', h, 'd');
        h(h <= 0) = nan;
        dat.v.lb.ld = sum(log(h(:)), 'omitnan');
        clear h

        % Lower bound
        % -----------
        % KL divergence between multivariate normal distributions
        % posterior: mean = r,    precision = H + (l+w)*L
        % prior:     mean = 0 ,   precision = (l+w)*L
        K = prod(opt.tpl.lat)*3;
        if opt.optimise.r.l
            loglam = spm_prob('Gamma', 'Elog', model.r.l, model.r.n, K);
        else
            loglam = log(model.r.l);
        end
        dat.v.lb.val = -0.5*( - K ...
                              - K * model.mixreg.w(2) * log(2*pi) ...
                              - K * model.mixreg.w(1) * loglam...
                              - model.mixreg.w(1) * opt.pg.ld ...
                              + dat.v.lb.ld ...
                              + model.mixreg.w(1) * model.r.l * dat.v.lb.reg ...
                              + model.mixreg.w(1) * model.r.l * dat.v.lb.tr );
        dat.v.lb.geod = -0.5*model.mixreg.w(2)*( ...
                                  K * log(2*pi) ...
                                - opt.pg.ld ...
                                + dat.v.lb.regv ...
                                + dat.v.lb.tr ...
                                + dat.v.lb.uncty );
        dat.v.lb.type = 'kl';
    else
        dat.v.lb.val   = 0;
        dat.v.lb.tr    = 0;
        dat.v.lb.uncty = 0;
        dat.v.lb.reg   = 0;
        dat.v.lb.regv  = 0;
        dat.v.lb.type  = '';
    end
        
    % ------------
    %    Latent
    % ------------
    if opt.optimise.z.z
        if opt.optimise.z.A
            logdetA = spm_prob('Wishart', 'ELogDet', model.z.A, model.z.n, 'normal');
        else
            logdetA = spm_matcomp('LogDet', model.z.A);
        end
        dat.z.lb.type = 'kl';
        dat.z.lb.val  = -0.5*( trace((dat.z.S + dat.z.zz) * model.z.A) ...
                               - logdetA ...
                               - spm_matcomp('LogDet', dat.z.S) ...
                               - opt.pg.K );
    else
        dat.z.lb.val  = 0;
        dat.z.lb.type = '';
    end
    
    % ------------
    %    Affine
    % ------------
    if opt.optimise.q.q && opt.q.Mr
        rind = opt.q.rind;
        qq   = dat.q.qq(rind,rind);
        Sq   = dat.q.S(rind,rind);
        dat.q.lb.type = 'kl';
        if opt.optimise.q.A
            logdetA = spm_prob('Wishart', 'ELogDet', model.q.A, model.q.n, 'normal');
        else
            logdetA = spm_matcomp('LogDet', model.q.A);
        end
        dat.q.lb.val  = -0.5*( trace((Sq + qq) * model.q.A) ...
                               - logdetA ...
                               + spm_matcomp('LogDet', Sq) ...
                               - opt.q.Mr );
    else
        dat.q.lb.type = '';
        dat.q.lb.val  = 0;
    end

    
    % --------------
    %    Mixture
    % --------------
    if opt.optimise.mixreg.w || opt.optimise.mixreg.a
        N = model.mixreg.n;
        A = model.mixreg.a;
        N0 = opt.mixreg.n0;
        A0 = opt.mixreg.a0;
        dg = spm_prob('DiGamma', N*A) - spm_prob('DiGamma', N);
        dig = spm_prob('DiGamma', N*(1-A)) - spm_prob('DiGamma', N);
    end
    if opt.optimise.mixreg.w
        W1 = model.mixreg.w(1);
        W2 = model.mixreg.w(2);
        model.lb.wr.val = W1*(log(W1) - dg) + W2*(log(W2) + dg);
        model.lb.wr.type = '-kl';
    end
    if opt.optimise.mixreg.a
        model.lb.ar.val =   (N*A - N0*A0) * dg ...
                          + (N*(1-A) - N0*(1-A0)) * dig...
                          + beta_norm(N0*A0,N0*(1-A0)) ...
                          - beta_norm(N*A,N*(1-A));
        model.lb.ar.type = '-kl';
    end

end

function b = beta_norm(a,b)
    b = gammaln(a) + gammaln(b) - gammaln(a+b);
end

function [dat, model] = batchInitLowerBound(dat, model, opt)

    % --- Detect parallelisation scheme
    if strcmpi(opt.split.loop, 'subject') && opt.split.par > 0
        batch = opt.split.batch;
    else
        batch = 1;
    end

    % Init log-likelihood
    % -------------------
    if isfield(model.lb, 'q'),   model.lb.q.val  = 0; end
    if isfield(model.lb, 'z'),   model.lb.z.val  = 0; end
    if isfield(model.lb, 'r'),   model.lb.r.val  = 0; end % KL Residual
    if isfield(model.lb, 'v'),   model.lb.v.val  = 0; end % LL "Geodesic"
    if isfield(model.r,  'tr'),  model.r.tr      = 0; end
    model.r.tr      = 0;
    model.r.uncty   = 0;
    model.r.reg     = 0;
    
    % --- Batch processing
    N = numel(dat);
    if opt.ui.verbose, before = plotBatchBegin('Init Lap'); end
    for i=1:ceil(N/batch)
        n1 = (i-1)*batch + 1;
        ne = min(N, i*batch);

        if opt.ui.verbose, before = plotBatch(i, batch, N, 50, before); end
    
        [opt.dist, dat(n1:ne)] = distribute(opt.dist, 'pgra_batch', 'OneInitLowerBound', 'inplace', dat(n1:ne), model, opt);
        
        for n=n1:ne
            
            % Add individual contributions
            % ----------------------------
            if isfield(model.lb, 'q'),    model.lb.q.val = model.lb.q.val + dat(n).q.lb.val;   end
            if isfield(model.lb, 'z'),    model.lb.z.val = model.lb.z.val + dat(n).z.lb.val;   end
            if isfield(model.lb, 'r'),    model.lb.r.val = model.lb.r.val + dat(n).v.lb.val;   end
            if isfield(model.lb, 'v'),    model.lb.v.val = model.lb.v.val + dat(n).v.lb.geod;  end
            if isfield(model.r,  'tr'),   model.r.tr     = model.r.tr     + dat(n).v.lb.tr;    end
            if isfield(model.r, 'uncty'), model.r.uncty  = model.r.uncty  + dat(n).v.lb.uncty; end
            if isfield(model.r, 'reg'),   model.r.reg    = model.r.reg    + dat(n).v.lb.reg;   end
        end
        
    end
    if opt.ui.verbose, plotBatchEnd; end
end

%% ------------------------------------------------------------------------
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
    
    % Init subjects
    % -------------
    nq = 0;
    if opt.ui.verbose, before = plotBatchBegin('Init Q'); end
    for n=1:numel(dat)
        nq = nq + 1;
        if opt.ui.verbose, before = plotBatch(nq, 1, opt.N, 50, before); end
        q           = create([opt.q.M, 1]);
        dat(n).q.q  = q;
        dat(n).q.qq = q*q';
        dat(n).q.S  = zeros(opt.q.M);
        dat(n).q.A  = exponentiateAffine(q, opt.q.B);
        model.q.q   = model.q.q  + dat(n).q.q;
        model.q.qq  = model.q.qq + dat(n).q.qq;
        dat(n).q.ok  = 1; % for GN failure tracking
        dat(n).q.ok2 = 0; % for GN failure tracking
    end
    if opt.ui.verbose, plotBatchEnd; end

end

%% ------------------------------------------------------------------------
%    Init Residual
% -------------------------------------------------------------------------

function dat = oneInitResidual(dat, model, opt, mode)

    % --- Detect parallelisation scheme
    if strcmpi(opt.split.loop, 'subject')
        loop = '';
        par  = 0;
    else
        loop = opt.split.loop;
        par  = opt.split.par;
    end
        
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
    dat.v.ok  = 1; % for GN failure tracking
    dat.v.ok2 = 0; % for GN failure tracking
    dat.v.v = prepareOnDisk(dat.v.v, [opt.tpl.lat 3]);
    dat.v.r = prepareOnDisk(dat.v.r, [opt.tpl.lat 3]);
    if strcmpi(mode, 'zero')
        dat.v.r(:)   = 0;
        dat.v.lb.reg = 0;
        r = 0;
    else
        r = create([opt.tpl.lat 3], 'single');
        dat.v.r(:,:,:,:) = r;
        m = spm_diffeo('vel2mom', r, double([opt.tpl.vs opt.pg.prm]));
        dat.v.lb.reg = r(:)' * m(:);
        clear m
    end
    if opt.pg.provided
        wz = reconstructVelocity('latent', dat.z.z, 'subspace', model.pg.w, ...
                                 'loop', loop, 'par', par);
        v = r + wz;
        clear wz
        dat.v.v(:,:,:,:) = v;
    else
        dat.v.v(:) = r;
        v = r;
    end
    clear r
    if ~opt.pg.provided && strcmpi(mode, 'zero')
        dat.v.lb.regv = 0;
    else
        m = spm_diffeo('vel2mom', v, double([opt.tpl.vs opt.pg.prm]));
        dat.v.lb.regv = v(:)'*m(:);
    end
    clear v

end

function [dat, model] = batchInitResidual(mode, dat, model, opt)

    % --- Detect parallelisation scheme
    if strcmpi(opt.split.loop, 'subject') && opt.split.par > 0
        batch = opt.split.batch;
    else
        batch = 1;
    end

    % Init log-likelihood
    % -------------------
    model.r.reg     = 0;
    
    % --- Batch processing
    N = numel(dat);
    if opt.ui.verbose, before = plotBatchBegin('Init Res'); end
    for i=1:ceil(N/batch)
        n1 = (i-1)*batch + 1;
        ne = min(N, i*batch);

        if opt.ui.verbose, before = plotBatch(i, batch, N, 50, before); end
    
        [opt.dist, dat(n1:ne)] = distribute(opt.dist, 'pgra_batch', 'OneInitResidual', 'inplace', dat(n1:ne), model, opt, mode);
        
    end
    if opt.ui.verbose, plotBatchEnd; end
end

%% ------------------------------------------------------------------------
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
    
    % Init subjects
    % -------------
    if opt.ui.verbose, before = plotBatchBegin('Init Z'); end
    for n=1:N
        if opt.ui.verbose, before = plotBatch(n, 1, N, 50, before); end
        z = create([opt.pg.K, 1]);
        dat(n).z.z   = z;
        dat(n).z.zz  = z*z';
        dat(n).z.S   = inv(model.mixreg.w(2) * model.pg.ww + model.z.A);
        dat(n).z.ok  = 1;
        dat(n).z.ok2 = 0;
        model.z.z    = model.z.z  + dat(n).z.z;
        model.z.zz   = model.z.zz + dat(n).z.zz;
        model.z.S    = model.z.S  + dat(n).z.S;
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
        end
        if opt.ui.verbose, plotBatchEnd; end
        
    end

end

%% ------------------------------------------------------------------------
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

%% ========================================================================
%    E-step
% =========================================================================

%% ------------------------------------------------------------------------
% Fit Affine
% -------------------------------------------------------------------------

function dat = oneFitAffine(dat, model, opt)

    % Detect parallelisation scheme
    % -----------------------------
    if strcmpi(opt.split.loop, 'subject')
        loop = '';
        par  = 0;
        verbose = false;
    else
        loop    = opt.split.loop;
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
            iphi = spm_warps('identity', opt.lat);
            phi  = [];
            jac  = [];
        end
    end
        
    % Penalise previous failure
    % -------------------------
    cumok = false;
    if opt.iter.pena && dat.q.ok < 0
        dat.q.ok = dat.q.ok + 1;
    else

        % Gauss-Newton iterations
        % -----------------------
        % It is useful to actually find a mode of the posterior (and not only
        % an improved value) when we use the Laplace precision for the update
        % of W. In that case, setting gnit > 1 might help converge faster.
        for i=1:opt.iter.gn

            % Compute gradient/hessian
            % ------------------------

            if strcmpi(opt.match, 'pull')
                [g, h] = ghMatchingAffine(noisemodel, ...
                    dat.tpl.wmu, dat.f.f, ...
                    model.tpl.gmu, dat.q.A, opt.q.B, ...
                    phi, jac, 'ipsi', dat.v.ipsi, ...
                    'Mmu', model.tpl.M, 'loop', loop, 'par', par, ...
                    'debug', opt.ui.debug, 'approx', opt.q.hapx);
            else
                [g, h] = ghMatchingAffine(noisemodel, ...
                    model.tpl.mu, dat.f.pf, ...
                    model.tpl.gmu, dat.q.A, opt.q.B, ...
                    phi, jac, 'count', dat.f.c, 'bb', dat.f.bb, ...
                    'Mmu', model.tpl.M, 'loop', loop, 'par', par, ...
                    'debug', opt.ui.debug, 'approx', opt.q.hapx);
            end

            if checkarray(model.q.A)
                rind = opt.q.rind;
                [gq, hq] = ghPriorAffine(dat.q.q(rind), model.q.A, 'debug', opt.debug);
                g(rind)      = g(rind)      + gq;
                h(rind,rind) = h(rind,rind) + hq;
                clear gq hq
            else
                rind = [];
            end

            h = spm_matcomp('LoadDiag', h); % Additional regularisation for robustness)

            % Compute search direction
            % ------------------------
            dq = -h\g;

            % Line search
            % -----------
            if opt.tpl.cat && strcmpi(opt.match, 'pull')
                a = model.tpl.a;
            else
                a = model.tpl.mu;
            end
            result = lsAffine(...
                noisemodel, dq, dat.q.q, dat.f.lb.val, a, dat.f.f, ...
                'B', opt.q.B, 'regq', model.q.A, 'rind', rind, ...
                'iphi', iphi, 'Mf', dat.f.M, 'Mmu', model.tpl.M, ...
                'match', opt.match, ...
                'nit', opt.iter.ls, 'loop', loop, 'par', par, ...
                'verbose', verbose, 'debug', opt.ui.debug, ...
                'pf', dat.f.pf, 'c', dat.f.c, 'wa', dat.tpl.wa, 'wmu', dat.tpl.wmu);
            
            % Store better values
            % -------------------
            cumok = cumok || result.ok;
            compute_hessian = result.ok;
            if result.ok
                dat.q.q       = result.q;
                dat.q.qq      = result.q * result.q';
                dat.q.A       = result.A;
                dat.f.lb.val  = result.llm;
                dat.v.ipsi    = copyarray(result.ipsi, dat.v.ipsi);
                dat.f.pf      = copyarray(result.pf,   dat.f.pf);
                dat.f.c       = copyarray(result.c,    dat.f.c);
                dat.f.bb      = result.bb;
                if strcmpi(opt.match, 'pull')
                    dat.tpl.wmu = copyarray(result.wmu, dat.tpl.wmu);
                    rmarray(result.wa);
                end
            else
                break
            end

        end % < GN iterations
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
        
    end % < penalise previous failure
    
    if cumok % < Only update LB if success
        
        if opt.q.Mr
            
            % Prior / KL-divergence
            % ---------------------

            if compute_hessian
                if strcmpi(opt.match, 'pull')
                    dat.q.h = ghMatchingAffine(noisemodel, ...
                        dat.tpl.wmu, dat.f.f, ...
                        model.tpl.gmu, dat.q.A, opt.q.B, ...
                        phi, jac, 'ipsi', dat.v.ipsi, ...
                        'hessian', true, ...
                        'Mmu', model.tpl.M, 'loop', loop, 'par', par, ...
                        'debug', opt.ui.debug, 'approx', opt.q.hapx);
                else
                    dat.q.h = ghMatchingAffine(noisemodel, ...
                        model.tpl.mu, dat.f.pf, ...
                        model.tpl.gmu, dat.q.A, opt.q.B, ...
                        phi, jac, 'count', dat.f.c, 'bb', dat.f.bb, ...
                        'hessian', true, ...
                        'Mmu', model.tpl.M, 'loop', loop, 'par', par, ...
                        'debug', opt.ui.debug, 'approx', opt.q.hapx);
                end

                if checkarray(model.q.A)
                    [~, hq] = ghPriorAffine(dat.q(rind), model.q.A, 'debug', opt.ui.debug);
                    dat.q.h(rind,rind) = dat.q.h(rind,rind) + hq;
                    clear hq
                end

                dat.q.h = spm_matcomp('LoadDiag', dat.q.h); % Additional regularisation for robustness
            end
            dat.q.S = inv(dat.q.h);
        else
            dat.q.S = 0;
        end
    
        % Lower bound
        % -----------
        if opt.q.Mr
            rind = opt.q.rind;
            qq   = dat.q.qq(rind,rind);
            Sq   = dat.q.S(rind,rind);
            dat.q.lb.type = 'kl';
            dat.q.lb.val  = -0.5*( trace((Sq + qq) * model.q.A) ...
                                   - spm_prob('Wishart', 'ELogDet', model.q.A, model.q.n, 'normal') ...
                                   + spm_matcomp('LogDet', Sq) ...
                                   - opt.q.Mr );
        end
                      
    end % < cumok
    
    % Cleaning (just in case)
    % --------
    dat.v.iphi = rmarray(dat.v.iphi);
    dat.v.phi  = rmarray(dat.v.phi);
    dat.v.jac  = rmarray(dat.v.jac);
    dat.tpl.wa = rmarray(dat.tpl.wa);
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
        [opt.dist, dat(n1:ne)] = distribute(opt.dist, 'pgra_batch', 'oneFitAffine', 'inplace', dat(n1:ne), model, opt);
        
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


%% ------------------------------------------------------------------------
% Fit Latent coordinates
% -------------------------------------------------------------------------

function dat = oneFitLatent(dat, model, opt)
    
    % Detect parallelisation scheme
    % -----------------------------
    if strcmpi(opt.split.loop, 'subject')
        loop = '';
        par  = 0;
        verbose = false;
    else
        loop    = opt.split.loop;
        par     = opt.split.par;
        verbose = opt.split.verbose;
    end
    
    if isfield(dat, 'model')
        noisemodel = dat.model;
    else
        noisemodel = opt.model;
    end
    
    if isfield(dat, 'q') && isfield(dat.q, 'A'),  A = dat.q.A;
    else,                                         A = eye(4);  end
    
    % Penalise previous failure
    % -------------------------        
    cumok = false;
    if opt.iter.pena && dat.z.ok < 0
        dat.z.ok = dat.z.ok + 1;
    else

        % Gauss-Newton iterations
        % -----------------------
        % It is useful to actually find a mode of the posterior (and not only
        % an improved value) when we use the Laplace precision for the update
        % of W. In that case, setting gnit > 1 might help converge faster.
        for i=1:opt.iter.gn

            % Compute gradient/hessian
            % ------------------------

            if strcmpi(opt.match, 'pull')
                [g, h] = ghMatchingLatent(noisemodel, ...
                    dat.tpl.wmu, dat.f.f, model.tpl.gmu, model.pg.w, ...
                    'ipsi', dat.v.ipsi, 'loop', loop, ...
                    'par', par, 'debug', opt.ui.debug);
            else
                [g, h] = ghMatchingLatent(noisemodel, ...
                    model.tpl.mu, dat.f.pf, model.tpl.gmu, model.pg.w, ...
                    'count', dat.f.c, 'bb', dat.f.bb, 'loop', loop, ...
                    'par', par, 'debug', opt.ui.debug);
            end

            [gz, hz] = ghPriorLatent(dat.z.z, model.z.A + model.mixreg.w(2) * model.pg.ww, 'debug', opt.ui.debug);
            g = g + gz;
            h = h + hz;
            clear gz hz

            % Part of geodesic prior
            if model.mixreg.w(1) > 0 && checkarray(dat.v.r)
                m = spm_diffeo('vel2mom', single(numeric(dat.v.r)), double([opt.tpl.vs opt.pg.prm]));
                for k=1:opt.pg.K
                    w1 = single(model.pg.w(:,:,:,:,k));
                    g(k) = g(k) + model.mixreg.w(1) * w1(:)' * m(:);
                end
                clear w1
            end
            clear m
            
            h = spm_matcomp('LoadDiag', h); % Additional regularisation for robustness

            % Compute search direction
            % ------------------------
            dz = -h \ g;

            % Line search
            % -----------
            if opt.tpl.cat && strcmpi(opt.match, 'pull')
                a = model.tpl.a;
            else
                a = model.tpl.mu;
            end
            result = lsLatent(...
                noisemodel, dz, dat.z.z, dat.v.v, dat.f.lb.val, ...
                model.pg.w, a, dat.f.f, ...
                'regz', model.z.A, 'geod', model.mixreg.w(2), ...
                'A', A, 'Mf', dat.f.M, 'Mmu', model.tpl.M, ...
                'nit', opt.iter.ls, 'itgr', opt.iter.itg, ...
                'prm', opt.pg.prm, 'bnd', opt.pg.bnd, ...
                'match', opt.match, ...
                'par', par, 'loop', loop, ...
                'verbose', verbose, 'debug', opt.ui.debug, ...
                'pf', dat.f.pf, 'c', dat.f.c, 'wa', dat.tpl.wa, 'wmu', dat.tpl.wmu);

            % Store better values
            % -------------------
            cumok = cumok || result.ok;
            compute_hessian = result.ok;
            if result.ok
                dat.z.z       = result.z;
                dat.z.zz      = result.z * result.z';
                dat.f.lb.val  = result.llm;
                dat.v.v       = copyarray(result.v,    dat.v.v);
                if model.mixreg.w(2)
                    m = spm_diffeo('vel2mom', result.v, double([opt.tpl.vs opt.pg.prm]));
                    dat.v.lb.regv = result.v(:)' * m(:);
                    clear m
                end
                dat.v.iphi    = copyarray(result.iphi, dat.v.iphi);
                dat.v.ipsi    = copyarray(result.ipsi, dat.v.ipsi);
                dat.f.pf      = copyarray(result.pf,   dat.f.pf);
                dat.f.c       = copyarray(result.c,    dat.f.c);
                dat.f.bb      = result.bb;
                if strcmpi(opt.match, 'pull')
                    dat.tpl.wmu = copyarray(result.wmu, dat.tpl.wmu);
                    rmarray(result.wa);
                end
            else
                break
            end
            clear result

        end % < GN iterations
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
        
    end % < penalise previous failure
    
    if cumok % < Only update LB if success
        
        % Prior / KL divergence
        % ---------------------

        % Hessian (Laplace approximation)
        % -------------------------------
        if compute_hessian
            if strcmpi(opt.match, 'pull')
                h = ghMatchingLatent(noisemodel, ...
                    dat.tpl.wmu, dat.f.f, model.tpl.gmu, model.pg.w, ...
                    'ipsi', dat.v.ipsi, 'hessian', true, ...
                    'loop', loop, 'par', par, 'debug', opt.ui.debug);
            else
                h = ghMatchingLatent(noisemodel, ...
                    model.tpl.mu, dat.f.pf, model.tpl.gmu, model.pg.w, ...
                    'count', dat.f.c, 'bb', dat.f.bb, 'hessian', true, ...
                    'loop', loop, 'par', par, 'debug', opt.ui.debug);
            end

            [~, hz] = ghPriorLatent(dat.z.z, model.z.A + model.mixreg.w(2) * model.pg.ww);
            h = h + hz;
            clear hz

            h = spm_matcomp('LoadDiag', h); % Additional regularisation for robustness

        end
        dat.z.S = inv(h);
        clear h
    
    end % cumok               
    
    % Cleaning
    % --------
    % Just in case
    dat.v.iphi = rmarray(dat.v.iphi);
    dat.tpl.wa = rmarray(dat.tpl.wa);
end

function [dat, model] = batchFitLatent(dat, model, opt)

    % --- Detect parallelisation scheme
    if strcmpi(opt.split.loop, 'subject') && opt.split.par > 0
        batch = opt.split.batch;
    else
        batch = 1;
    end

    % Init gradient/hessian
    % ---------------------
    model.z.z      = zeros(opt.pg.K, 1);
    model.z.zz     = zeros(opt.pg.K);
    model.z.S      = zeros(opt.pg.K);
    
    N   = numel(dat);
    Nok = 0;
    for n=1:N
        if isfield(dat(n).z, 'ok')
            Nok = Nok + (dat(n).z.ok >= 0);
        end
    end
    
    % --- Batch processing
    if opt.ui.verbose, before = plotBatchBegin('Fit Latent'); end
    for i=1:ceil(N/batch)
        n1 = (i-1)*batch + 1;
        ne = min(N, i*batch);

        if opt.ui.verbose, before = plotBatch(i, batch, N, 50, before); end
    
        % Compute subjects grad/hess w.r.t. initial velocity
        % --------------------------------------------------
        [opt.dist, dat(n1:ne)] = distribute(opt.dist, 'pgra_batch', 'oneFitLatent', 'inplace', dat(n1:ne), model, opt);
        
        for n=n1:ne
            
            % Add individual contributions
            % ----------------------------
            model.z.z      = model.z.z      + dat(n).z.z;
            model.z.zz     = model.z.zz     + dat(n).z.zz;
            model.z.S      = model.z.S      + dat(n).z.S;
        end
        
    end
    okz = 0;
    for n=1:N
        if isfield(dat(n).z, 'ok')
            okz = okz + (dat(n).z.ok > 0);
        end
    end
    if opt.ui.verbose, fprintf(' | %3d / %3d / %3d', okz, Nok, N); end
    if opt.ui.verbose, plotBatchEnd; end
end


%% ------------------------------------------------------------------------
% Fit Residual field
% -------------------------------------------------------------------------

function dat = oneFitResidual(dat, model, opt)

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
    cumok = false;
    if opt.iter.pena && dat.v.ok < 0
        dat.v.ok = dat.v.ok + 1;
    else
        
        % Gauss-Newton iterations
        % -----------------------
        for i=1:opt.iter.gn
            
            % Compute gradient/hessian
            % ------------------------

            if strcmpi(opt.match, 'pull')
                [g, h] = ghMatchingVel(...
                    noisemodel, ...
                    dat.tpl.wmu, dat.f.f, model.tpl.gmu, 'ipsi', dat.v.ipsi, ...
                    'loop', loop, 'par', par, 'debug', opt.ui.debug);
            else
                g = zeros([opt.tpl.lat 3], 'single');
                h = zeros([opt.tpl.lat 6], 'single');
                [g(dat.f.bb.x,dat.f.bb.y,dat.f.bb.z,:), ...
                 h(dat.f.bb.x,dat.f.bb.y,dat.f.bb.z,:)] ...
                    = ghMatchingVel(noisemodel, ...
                        model.tpl.mu, dat.f.pf, model.tpl.gmu, ...
                        'count', dat.f.c, 'bb', dat.f.bb, ...
                        'loop', loop, 'par', par, 'debug', opt.ui.debug);
            end
            v = numeric(dat.v.v);
            r = numeric(dat.v.r);
            g = g + ghPriorVel(r, opt.tpl.vs, (model.mixreg.w(1) * model.r.l + model.mixreg.w(2)) * opt.pg.prm, opt.pg.bnd);

            % Part of geodesic prior
            if model.mixreg.w(2)
               g = g + spm_diffeo('vel2mom', single(v-r), double([opt.tpl.vs model.mixreg.w(2) * opt.pg.prm]));
            end
            
            % Compute search direction
            % ------------------------
            dv = -spm_diffeo('fmg', single(h), single(g), ...
                double([opt.tpl.vs (model.mixreg.w(1) * model.r.l + model.mixreg.w(2)) * opt.pg.prm 2 2]));
            clear g

            % Line search
            % -----------
            if opt.tpl.cat && strcmpi(opt.match, 'pull')
                a = model.tpl.a;
            else
                a = model.tpl.mu;
            end
            result = lsVelocity(...
                noisemodel, dv, r, dat.f.lb.val, a, dat.f.f, ...
                'v0', v, 'lam', model.mixreg.w(1) * model.r.l, 'geod', model.mixreg.w(2), ...
                'prm', opt.pg.prm, 'itgr', opt.iter.itg, 'bnd', opt.pg.bnd, ...
                'A', A, 'Mf', dat.f.M, 'Mmu', model.tpl.M, ...
                'match', opt.match, ...
                'nit', opt.iter.ls,  'par', par, 'loop', loop, ...
                'verbose', verbose, 'debug', opt.ui.debug, ...
                'pf', dat.f.pf, 'c', dat.f.c, 'wa', dat.tpl.wa, 'wmu', dat.tpl.wmu);

            % Store better values
            % -------------------
            cumok = cumok || result.ok;
            compute_hessian = result.ok;
            if result.ok
                dat.f.lb.val  = result.match;
                dat.v.v       = copyarray(result.v,    dat.v.v);
                if model.mixreg.w(2)
                    m = spm_diffeo('vel2mom', result.v, double([opt.tpl.vs opt.pg.prm]));
                    dat.v.lb.regv = result.v(:)' * m(:);
                    clear m
                end
                dat.v.r       = copyarray(result.r,    dat.v.r);
                dat.v.ipsi    = copyarray(result.ipsi, dat.v.ipsi);
                dat.f.pf      = copyarray(result.pf,   dat.f.pf);
                dat.f.c       = copyarray(result.c,    dat.f.c);
                dat.f.bb      = result.bb;
                if strcmpi(opt.match, 'pull')
                    dat.tpl.wmu   = copyarray(result.wmu, dat.tpl.wmu);
                    rmarray(result.wa);
                end
                r    = result.r;
                ipsi = result.ipsi;
            else
                break
            end
            clear result

        end % < GN iterations
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
    
    end % < penalise previous failure
    
    if cumok % < Only update LB if success
        
        % -----------
        % Lower bound
        % -----------
        
        % Regularisation part
        % -------------------
        m = spm_diffeo('vel2mom', r, double([opt.tpl.vs opt.pg.prm]));
        dat.v.lb.reg = r(:)' * m(:);
        clear r m
        
        % Trace(P\L) part
        % ---------------
        if compute_hessian
            if strcmpi(opt.match, 'pull')
                h = ghMatchingVel(noisemodel, ...
                    dat.tpl.wmu, dat.f.f, model.tpl.gmu, 'ipsi', ipsi, ...
                    'hessian', true, 'loop', loop, 'par', par, ...
                    'debug', opt.ui.debug);
            else
                h = zeros([opt.tpl.lat 6], 'single');
                h(dat.f.bb.x,dat.f.bb.y,dat.f.bb.z,:) = ghMatchingVel(...
                    noisemodel, ...
                    model.tpl.mu, dat.f.pf, model.tpl.gmu, ...
                    'count', dat.f.c, 'bb', dat.f.bb, ...
                    'hessian', true, 'loop', loop, 'par', par, ...
                    'debug', opt.ui.debug);
            end
        end
        clear ipsi
        dat.v.lb.tr = spm_diffeo('trapprox', h, double([opt.tpl.vs (model.mixreg.w(1)*model.r.l + model.mixreg.w(2)) * opt.pg.prm]));
        dat.v.lb.tr = dat.v.lb.tr(1);
        dat.v.lb.tr = dat.v.lb.tr / (model.mixreg.w(1)*model.r.l + model.mixreg.w(2));
        
        % LogDet(P)
        % ---------
        % Approximation where all off-diagonal elements of L are zero
        h(:,:,:,1) = h(:,:,:,1) * (model.mixreg.w(1)*model.r.l + model.mixreg.w(2)) * opt.pg.ker(1);
        h(:,:,:,2) = h(:,:,:,2) * (model.mixreg.w(1)*model.r.l + model.mixreg.w(2)) * opt.pg.ker(2);
        h(:,:,:,3) = h(:,:,:,3) * (model.mixreg.w(1)*model.r.l + model.mixreg.w(2)) * opt.pg.ker(3);
        if opt.model.dim == 2
            h(:,:,:,3) = 1;
        end
        h = spm_matcomp('Pointwise3', h, 'd');
        h(h <= 0) = nan;
        dat.v.lb.ld = sum(log(h(:)), 'omitnan');
        clear h
        
        % KL divergence
        % -------------
        K = prod(opt.tpl.lat)*3;
        if opt.optimise.r.l
            loglam = spm_prob('Gamma', 'Elog', model.r.l, model.r.n, K);
        else
            loglam = log(model.r.l);
        end
        dat.v.lb.val = -0.5*( - K ...
                              - K * model.mixreg.w(2) * log(2*pi) ...
                              - K * model.mixreg.w(1) * loglam...
                              - model.mixreg.w(1) * opt.pg.ld ...
                              + dat.v.lb.ld ...
                              + model.mixreg.w(1) * model.r.l * dat.v.lb.reg ...
                              + model.mixreg.w(1) * model.r.l * dat.v.lb.tr );
        dat.v.lb.geod = -0.5*model.mixreg.w(2)*( ...
                                  K * log(2*pi) ...
                                - opt.pg.ld ...
                                + dat.v.lb.regv ...
                                + dat.v.lb.tr ...
                                + dat.v.lb.uncty );
    
    end % < cumok
    
    % Cleaning
    % --------
    % Just in case
    dat.v.iphi = rmarray(dat.v.iphi);
    dat.tpl.wa = rmarray(dat.tpl.wa);
end

function [dat, model] = batchFitResidual(dat, model, opt)

    % --- Detect parallelisation scheme
    if strcmpi(opt.split.loop, 'subject') && opt.split.par > 0
        batch = opt.split.batch;
    else
        batch = 1;
    end

    % Init gradient/hessian
    % ---------------------
    model.lb.m.val  = 0;
    model.lb.r.val  = 0;
    model.lb.v.val  = 0;
    model.r.tr      = 0;
    model.r.reg     = 0;
    
    total = opt.N;
    N     = numel(dat);
    okpre = 0;
    for n=1:N
        if isfield(dat(n).v, 'ok')
            okpre = okpre + (dat(n).v.ok >= 0);
        end
    end
    
    % --- Batch processing
    if opt.ui.verbose, before = plotBatchBegin('Fit Res'); end
    for i=1:ceil(N/batch)
        n1 = (i-1)*batch + 1;
        ne = min(N, i*batch);

        if opt.ui.verbose, before = plotBatch(i, batch, N, 50, before); end
    
        % Compute subjects grad/hess w.r.t. initial velocity
        % --------------------------------------------------
        [opt.dist, dat(n1:ne)] = distribute(opt.dist, 'pgra_batch', 'oneFitResidual', 'inplace', dat(n1:ne), model, opt);
        
        
        for n=n1:ne
            
            % Add individual contributions
            % ----------------------------
            model.r.reg     = model.r.reg    + dat(n).v.lb.reg;
            model.r.tr      = model.r.tr     + dat(n).v.lb.tr;
            model.lb.r.val  = model.lb.r.val + dat(n).v.lb.val;
            model.lb.m.val  = model.lb.m.val + dat(n).f.lb.val;
            model.lb.v.val  = model.lb.v.val + dat(n).v.lb.geod;
            
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

%% ------------------------------------------------------------------------
%    M-step
% -------------------------------------------------------------------------

function dat = oneGradHessVelocity(dat, model, opt)

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
    
    if strcmpi(opt.match, 'pull')
        [g, h] = ghMatchingVel(noisemodel, ...
            dat.tpl.wmu, dat.f.f, model.tpl.gmu, 'ipsi', dat.v.ipsi, ...
            'loop', loop, 'par', par, 'debug', opt.ui.debug);
    else
        g = zeros([opt.tpl.lat 3], 'single');
        h = zeros([opt.tpl.lat 6], 'single');
        [g(dat.f.bb.x,dat.f.bb.y,dat.f.bb.z,:), ...
         h(dat.f.bb.x,dat.f.bb.y,dat.f.bb.z,:)] ...
            = ghMatchingVel(noisemodel, ...
                model.tpl.mu, dat.f.pf, model.tpl.gmu, ...
                'count', dat.f.c, 'bb', dat.f.bb, ...
                'loop', loop, 'par', par, 'debug', opt.ui.debug);
    end
    dat.v.g = copyarray(g, dat.v.g);
    dat.v.h = copyarray(h, dat.v.h);
end

% --------
% Subspace
% --------

function [dat, model] = batchGradHessSubspace(dat, model, opt)

    % --- Detect parallelisation scheme
    if strcmpi(opt.split.loop, 'subject') && opt.split.par > 0
        batch = opt.split.batch;
    else
        batch = 1;
    end

    % Init gradient/hessian
    % ---------------------
    dim        = [size(model.pg.w) 1 1 1];
    model.pg.g = prepareOnDisk(model.pg.g, [dim(1:3) 3 dim(5)], 'type', 'float32');
    model.pg.h = prepareOnDisk(model.pg.h, [dim(1:3) 6 dim(5)], 'type', 'float32');
    for k=1:dim(5)
        model.pg.g(:,:,:,:,k) = 0;
        model.pg.h(:,:,:,:,k) = 0;
    end
    
    % --- Batch processing
    if opt.ui.verbose, before = plotBatchBegin('GH PG'); end
    for i=1:ceil(opt.N/batch)
        n1 = (i-1)*batch + 1;
        ne = min(opt.N, i*batch);

        if opt.ui.verbose, before = plotBatch(i, batch, opt.N, 50, before); end
    
        % Compute subjects grad/hess w.r.t. initial velocity
        % --------------------------------------------------
        [opt.dist, dat(n1:ne)] = distribute(opt.dist, 'pgra_batch', 'oneGradHessVelocity', 'inplace', dat(n1:ne), model, opt);
        
        
        % Add individual contributions
        % ----------------------------
        if model.mixreg.w(2)
            for k=1:opt.pg.K
                gw = model.pg.g(:,:,:,:,k);
                hw = model.pg.h(:,:,:,:,k);
                for n=n1:ne
                    gv = numeric(dat(n).v.g);
                    hv = numeric(dat(n).v.h);
                    gw = gw + gv * model.mixreg.w(2) * single(dat(n).z.z(k));
                    hw = hw + hv * model.mixreg.w(2) * single(dat(n).z.z(k))^2;
                    clear gv hv
                end
                model.pg.g(:,:,:,:,k) = gw;
                model.pg.h(:,:,:,:,k) = hw;
                clear gw hw
            end
        end
        
        % Clear individual grad/hess
        % --------------------------
        for n=n1:ne
            dat(n).v.g = rmarray(dat(n).v.g);
            dat(n).v.h = rmarray(dat(n).v.h);
        end
        
    end
    if opt.ui.verbose, plotBatchEnd; end

    % Total residual momentum
    % -----------------------
    if model.mixreg.w(2)
        if opt.ui.verbose, before = plotBatchBegin('GH PGr'); end
        m = zeros([opt.tpl.lat 3], 'single');
        for n=1:opt.N
            if opt.ui.verbose, before = plotBatch(n, 1, opt.N, 50, before); end
            m = m + single(numeric(dat(n).v.r));
        end
        m = spm_diffeo('vel2mom', m, double([opt.tpl.vs opt.pg.prm]));
        if opt.ui.verbose, plotBatchEnd; end
    else
        m = 0;
    end
                        
    % Regularisation gradient
    % -----------------------
    reg = model.mixreg.w(2) * (model.z.zz + model.z.S) + opt.N * eye(opt.pg.K);
    for k=1:opt.pg.K
        lw = spm_diffeo('vel2mom', single(model.pg.w(:,:,:,:,k)), [opt.tpl.vs, opt.pg.prm]);
        model.pg.g(:,:,:,:,k) = model.pg.g(:,:,:,:,k) + reg(k,k) * lw + model.mixreg.w(2) * model.z.z(k) * m;
        clear lw
    end
    clear m
    
end

%% ------------------------------------------------------------------------
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
        
    switch lower(which)
        
        case 'matching'
            if strcmpi(opt.match, 'pull')
                dat.f.lb.val = llMatching(noisemodel, dat.tpl.wmu, dat.f.f, ...
                    'par', par, 'loop', loop, 'debug', opt.ui.debug);
            else
                dat.f.lb.val = llMatching(noisemodel, model.tpl.mu, dat.f.pf, dat.f.c, ...
                    'bb', dat.f.bb, 'par', par, 'loop', loop, 'debug', opt.ui.debug);
            end
            
        case 'template'
            if strcmpi(opt.match, 'pull')
                if opt.tpl.cat
                    dat.tpl.wa = pullTemplate(dat.v.ipsi, model.tpl.a, ...
                        'par', par, 'output', dat.tpl.wa, 'debug', opt.ui.debug);
                    dat.tpl.wmu = reconstructProbaTemplate(dat.tpl.wa, ...
                        'loop', loop, 'par', par, 'output', dat.tpl.wmu, ...
                        'debug', opt.ui.debug);
                    dat.tpl.wa = rmarray(dat.tpl.wa);
                else
                    dat.tpl.wmu = pullTemplate(dat.v.ipsi, model.tpl.mu, ...
                        'par', par, 'output', dat.tpl.wmu, 'debug', opt.ui.debug);
                end
                dat.f.lb.val = llMatching(noisemodel, dat.tpl.wmu, dat.f.f, ...
                    'par', par, 'loop', loop, 'debug', opt.ui.debug);
            else
                dat.f.lb.val = llMatching(noisemodel, model.tpl.mu, dat.f.pf, dat.f.c, ...
                    'bb', dat.f.bb, 'par', par, 'loop', loop, 'debug', opt.ui.debug);
            end
            
        case 'precisionz'
            if opt.optimise.z.z
                if opt.optimise.z.A
                    logdetA = spm_prob('Wishart', 'ELogDet', model.z.A, model.z.n, 'normal');
                else
                    logdetA = spm_matcomp('LogDet', model.z.A);
                end
                dat.z.lb.val  = -0.5*( trace((dat.z.S + dat.z.zz) * model.z.A) ...
                                       - logdetA ...
                                       - spm_matcomp('LogDet', dat.z.S) ...
                                       - opt.pg.K );
            end
    
        case 'precisionq'
            if opt.f.observed && opt.optimise.q.q && opt.q.Mr
                rind = opt.q.rind;
                qq   = dat.q.qq(rind,rind);
                Sq   = dat.q.S(rind,rind);
                dat.q.lb.type = 'kl';
                if opt.optimise.q.A
                    logdetA = spm_prob('Wishart', 'ELogDet', model.q.A, model.q.n, 'normal');
                else
                    logdetA = spm_matcomp('LogDet', model.q.A);
                end
                dat.q.lb.val  = -0.5*( trace((Sq + qq) * model.q.A) ...
                                       - logdetA ...
                                       + spm_matcomp('LogDet', Sq) ...
                                       - opt.q.Mr );
            end
                              
        case 'orthogonalise'
            if opt.optimise.z.z
                if opt.optimise.z.A
                    logdetA = spm_prob('Wishart', 'ELogDet', model.z.A, model.z.n, 'normal');
                else
                    logdetA = spm_matcomp('LogDet', model.z.A);
                end
                dat.z.lb.val  = -0.5*( trace((dat.z.S + dat.z.zz) * model.z.A) ...
                                       - logdetA ...
                                       - spm_matcomp('LogDet', dat.z.S) ...
                                       - opt.pg.K );
            end
                               
        case 'subspace'
            % I only need to reconstruct velocity of line search failed.
            % I should find a way to avoid doing it when not necessary
            v = reconstructVelocity('latent', dat.z.z, ...
                                    'subspace', model.pg.w, ...
                                    'residual', dat.v.r, ...
                                    'loop', loop, 'par', par);
            iphi = exponentiateVelocity(v, 'iphi', ....
                                        'itgr', opt.iter.itg, 'vs', opt.tpl.vs, ...
                                        'prm', opt.pg.prm, 'bnd', opt.pg.bnd);
            dat.v.v = copyarray(v, dat.v.v);
            m = spm_diffeo('vel2mom', v, double([opt.tpl.vs opt.pg.prm]));
            dat.v.lb.regv = v(:)'*m(:);
            clear m v
            if isfield(dat, 'q') && isfield(dat.q, 'A'),  A = dat.q.A;
            else,                                         A = eye(4);  end
            ipsi = reconstructIPsi(A, iphi, 'lat', opt.tpl.lat, ...
                                   'Mf',  dat.f.M, 'Mmu', model.tpl.M, ...
                                   'debug', opt.ui.debug);
            clear iphi
            [dat.f.pf, dat.f.c, dat.f.bb] = pushImage(ipsi, dat.f.f, opt.tpl.lat, ...
                                                      'output', {dat.f.pf, dat.f.c}, ...
                                                      'loop', loop, 'par', par, ...
                                                      'debug', opt.ui.debug);
            if strcmpi(opt.match, 'pull')
                if opt.tpl.cat
                    dat.tpl.wa = pullTemplate(ipsi, model.tpl.a, ...
                        'par', par, 'output', dat.tpl.wa, 'debug', opt.ui.debug);
                    dat.tpl.wmu = reconstructProbaTemplate(dat.tpl.wa, ...
                        'loop', loop, 'par', par, 'output', dat.tpl.wmu, ...
                        'debug', opt.ui.debug);
                    dat.tpl.wa = rmarray(dat.tpl.wa);
                else
                    dat.tpl.wmu = pullTemplate(dat.v.ipsi, model.tpl.mu, ...
                        'par', par, 'output', dat.tpl.wmu, 'debug', opt.ui.debug);
                end
                dat.v.ipsi = copyarray(ipsi, dat.v.ipsi);
                clear ipsi
                dat.f.lb.val = llMatching(noisemodel, dat.tpl.wmu, dat.f.f, ...
                    'par', par, 'loop', loop, 'debug', opt.ui.debug);
            else
                dat.f.lb.val = llMatching(noisemodel, model.tpl.mu, dat.f.pf, dat.f.c, ...
                    'bb', dat.f.bb, 'par', par, 'loop', loop, 'debug', opt.ui.debug);
            end
            K = prod(opt.tpl.lat) * 3;
            dat.v.lb.uncty = trace(dat.z.S * model.pg.ww);
            if opt.optimise.z.A
                logdetA = spm_prob('Wishart', 'ELogDet', model.z.A, model.z.n, 'normal');
            else
                logdetA = spm_matcomp('LogDet', model.z.A);
            end
            dat.z.lb.val  = -0.5*( trace((dat.z.S + dat.z.zz) * model.z.A) ...
                                   - logdetA ...
                                   - spm_matcomp('LogDet', dat.z.S) ...
                                   - opt.pg.K );
            if opt.optimise.r.l
                loglam = spm_prob('Gamma', 'Elog', model.r.l, model.r.n, K);
            else
                loglam = log(model.r.l);
            end
            dat.v.lb.val = -0.5*( - K ...
                                  - K * model.mixreg.w(2) * log(2*pi) ...
                                  - K * model.mixreg.w(1) * loglam...
                                  - model.mixreg.w(1) * opt.pg.ld ...
                                  + dat.v.lb.ld ...
                                  + model.mixreg.w(1) * model.r.l * dat.v.lb.reg ...
                                  + model.mixreg.w(1) * model.r.l * dat.v.lb.tr );
            dat.v.lb.geod = -0.5*model.mixreg.w(2)*( ...
                                      K * log(2*pi) ...
                                    - opt.pg.ld ...
                                    + dat.v.lb.regv ...
                                    + dat.v.lb.tr ...
                                    + dat.v.lb.uncty );
                               
        case 'latent'
            K = prod(opt.tpl.lat) * 3;
            dat.v.lb.uncty = trace(dat.z.S * model.pg.ww);
            if opt.optimise.z.A
                logdetA = spm_prob('Wishart', 'ELogDet', model.z.A, model.z.n, 'normal');
            else
                logdetA = spm_matcomp('LogDet', model.z.A);
            end
            dat.z.lb.val  = -0.5*( trace((dat.z.S + dat.z.zz) * model.z.A) ...
                                   - logdetA ...
                                   - spm_matcomp('LogDet', dat.z.S) ...
                                   - opt.pg.K );
            dat.v.lb.geod = -0.5*model.mixreg.w(2)*( ...
                                      K * log(2*pi) ...
                                    - opt.pg.ld ...
                                    + dat.v.lb.regv ...
                                    + dat.v.lb.tr ...
                                    + dat.v.lb.uncty );
            
        case 'lambda'
            K = prod(opt.tpl.lat) * 3;
            if opt.optimise.r.l
                loglam = spm_prob('Gamma', 'Elog', model.r.l, model.r.n, K);
            else
                loglam = log(model.r.l);
            end
            dat.v.lb.val = -0.5*( - K ...
                                  - K * model.mixreg.w(2) * log(2*pi) ...
                                  - K * model.mixreg.w(1) * loglam...
                                  - model.mixreg.w(1) * opt.pg.ld ...
                                  + dat.v.lb.ld ...
                                  + model.mixreg.w(1) * model.r.l * dat.v.lb.reg ...
                                  + model.mixreg.w(1) * model.r.l * dat.v.lb.tr );

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
    if isfield(model.lb, 'm'),   model.lb.m.val  = 0; end
    if isfield(model.lb, 'z'),   model.lb.z.val  = 0; end
    if isfield(model.lb, 'q'),   model.lb.q.val  = 0; end
    if isfield(model.lb, 'r'),   model.lb.r.val  = 0; end
    if isfield(model.lb, 'v'),   model.lb.v.val  = 0; end
    if isfield(model.lb, 'g'),   model.lb.g.val  = 0; end
    if isfield(model.r,  'tr'),  model.r.tr      = 0; end
    if isfield(model.r,  'reg'), model.r.reg     = 0; end
    
    N = numel(dat);
    
    % --- Batch processing
    if opt.ui.verbose, before = plotBatchBegin(['LB ' which(1:min(numel(which),3))]); end
    for i=1:ceil(N/batch)
        n1 = (i-1)*batch + 1;
        ne = min(N, i*batch);

        if opt.ui.verbose, before = plotBatch(i, batch, N, 50, before); end
    
        % Compute subjects grad/hess w.r.t. initial velocity
        % --------------------------------------------------
        [opt.dist, dat(n1:ne)] = distribute(opt.dist, 'pgra_batch', 'oneLB', 'inplace', dat(n1:ne), model, opt, which);
        
        for n=n1:ne
            
            % Add individual contributions
            % ----------------------------
            if isfield(model.lb, 'm'),   model.lb.m.val  = model.lb.m.val + dat(n).f.lb.val; end
            if isfield(model.lb, 'z'),   model.lb.z.val  = model.lb.z.val + dat(n).z.lb.val; end
            if isfield(model.lb, 'q'),   model.lb.q.val  = model.lb.q.val + dat(n).q.lb.val; end
            if isfield(model.lb, 'r'),   model.lb.r.val  = model.lb.r.val + dat(n).v.lb.val; end
            if isfield(model.lb, 'v'),   model.lb.v.val  = model.lb.v.val + dat(n).v.lb.geod; end
            if isfield(model.lb, 'g'),   model.lb.g.val  = model.lb.g.val + dat(n).v.lb.geod; end
            if isfield(model.r,  'tr'),  model.r.tr      = model.r.tr     + dat(n).v.lb.tr;  end
            if isfield(model.r,  'reg'), model.r.reg     = model.r.reg    + dat(n).v.lb.reg; end
            
        end
        
    end
    if opt.ui.verbose, plotBatchEnd; end

    % Model specific parts
    % --------------------
    switch lower(which)
        case 'orthogonalise'
            if isfield(model.lb, 'w')
                model.lb.w.val = llPriorSubspace(model.pg.w, opt.N * model.pg.ww, opt.pg.ld + prod(opt.tpl.lat)*3*log(opt.N));
            end
            if opt.z.n0 && isfield(model.lb, 'Az')
                model.lb.Az.val = -spm_prob('Wishart', 'kl', ...
                                           model.z.A,   model.z.n, ...
                                           opt.z.A0,    opt.z.n0, ...
                                           'normal');
            end
        case 'subspace'
            if isfield(model.lb, 'w')
                model.lb.w.val = llPriorSubspace(model.pg.w, model.pg.n * model.pg.ww, opt.pg.ld + prod(opt.tpl.lat)*3*log(model.pg.n));
            end
        case 'precisionz'
            if opt.z.n0 && isfield(model.lb, 'Az')
                model.lb.Az.val = -spm_prob('Wishart', 'kl', ...
                                            model.z.A,   model.z.n, ...
                                            opt.z.A0,    opt.z.n0, ...
                                            'normal');
            end
        case 'precisionq'
            if opt.q.n0 && opt.q.Mr && isfield(model.lb, 'Aq')
                model.lb.Aq.val = -spm_prob('Wishart', 'kl', ...
                                            model.q.A,   model.q.n, ...
                                            opt.q.A0,    opt.q.n0, ...
                                            'normal');
            end
        case 'lambda'
            if opt.r.n0 && isfield(model.lb, 'l')
                model.lb.l.val  = -spm_prob('Gamma', 'kl', ...
                                            model.r.l, model.r.n, ...
                                            opt.r.l0,  opt.r.n0, ...
                                            prod(opt.tpl.lat)*3, 'normal');
            end
    end
end
