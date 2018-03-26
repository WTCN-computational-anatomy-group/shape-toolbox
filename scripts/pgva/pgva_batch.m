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
% FORMAT [dat, model] = pgva_batch('InitLaplace', dat, model, opt)
%   Initialise covariance for Laplace approximations + lower bound
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
% FORMAT dat = pgva_batch('OneInitLaplace',  dat, model, opt)
% FORMAT dat = pgva_batch('OneFitAffine',    dat, model, opt)
% FORMAT dat = pgva_batch('OneFitLatent',    dat, model, opt)
% FORMAT dat = pgva_batch('OneFitVelocity',  dat, model, opt)
% FORMAT dat = pgva_batch('OneLB',           dat, model, opt, var)
% FORMAT dat = pgva_batch('OneMomentum',     dat, model, opt)
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
        case 'initlowerbound'
            [varargout{1:nargout}] = batchInitLowerBound(varargin{:});
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
        case 'oneinitlowerbound'
            [varargout{1:nargout}] = oneInitLowerBound(varargin{:});
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
%    Momentum
% -------------------------------------------------------------------------

function dat = oneMomentum(dat, opt)

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
    
        [opt.dist, dat(n1:ne)] = distribute(opt.dist, 'pgva_batch', 'OneMomentum', 'inplace', dat(n1:ne), opt);
        
    end
    if opt.ui.verbose, plotBatchEnd; end
end

%% ------------------------------------------------------------------------
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
    
        [opt.dist, dat(n1:ne)] = distribute(opt.dist, 'pgva_batch', 'OneInitPush', 'inplace', dat(n1:ne), model, opt);
        
    end
    if opt.ui.verbose, plotBatchEnd; end
end

%% ------------------------------------------------------------------------
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
    if opt.tpl.cat
        dat.tpl.scale = ones(1, opt.model.nc);
        dat.f.ss0 = 0;
        dat.f.ss1 = zeros(1, opt.model.nc);
        for z=1:size(dat.f.f, 3)
            f1        = dat.f.f(:,:,z,:);
            msk       = all(isfinite(f1), 4);
            dat.f.ss0 = dat.f.ss0 + sumall(msk);
            for k=1:opt.model.nc
                fk = f1(:,:,:,k);
                dat.f.ss1(k) = dat.f.ss1(k) + sum(fk(msk));
            end
            clear fk f1 msk
        end
        dat.tpl.wa = pullTemplate(dat.v.ipsi, model.tpl.a, ...
            'par', par, 'output', dat.tpl.wa, 'debug', opt.ui.debug);
        [dat.tpl.wmu, ssmu] = reconstructProbaTemplate(dat.tpl.wa, ...
            'scale', dat.tpl.scale, ...
            'loop', loop, 'par', par, 'output', dat.tpl.wmu, ...
            'debug', opt.ui.debug);
        if opt.optimise.tpl.scale
            % Individual probability scaling
            % ------------------------------
            for i=1:opt.iter.scl
                dimf = [size(dat.f.f) 1];
                dat.tpl.scale = (dat.f.ss1(:) * prod(dimf(1:3)) ./ (ssmu(:) * dat.f.ss0));
                [dat.tpl.wmu, ssmu] = reconstructProbaTemplate(dat.tpl.wa, ...
                    'scale', dat.tpl.scale, ...
                    'loop', loop, 'par', par, 'output', dat.tpl.wmu, ...
                    'debug', opt.ui.debug);
            end
        end
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
    
        [opt.dist, dat(n1:ne)] = distribute(opt.dist, 'pgva_batch', 'OneInitPull', 'inplace', dat(n1:ne), model, opt);
        
        for n=n1:ne
            if dat(n).f.observed
                model.lb.m.val = model.lb.m.val + dat(n).f.lb.val;
            end
        end
    end
    if opt.ui.verbose, plotBatchEnd; end
end
%% ------------------------------------------------------------------------
%    Init Lower Bound
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
    %    Velocity
    % --------------
    if dat.v.observed
        
        % Lower bound
        % -----------
        % log-likelihood of a mixture of two Gaussian distributions
        K = prod(opt.tpl.lat)*3;
        if opt.optimise.v.l
            loglam = spm_prob('Gamma', 'Elog', model.v.l, model.v.n, K);
        else
            loglam = log(model.v.l);
        end
        dat.v.lb.uncty = trace(dat.z.S * model.pg.ww);
        dat.v.lb.ll1 = -0.5*( K*log(2*pi) ...
                              - K * sloglam ...
                              - opt.pg.ld ...
                              + model.v.l * dat.v.lb.reg ...
                              + model.v.l * dat.v.lb.uncty );
        dat.v.lb.ll2 = -0.5*( K * log(2*pi) ...
                              - opt.pg.ld ...
                              + dat.v.lb.regv );
        dat.v.lb.val = model.mixreg.w(1) * dat.v.lb.ll1 + ...
                       model.mixreg.w(2) * dat.v.lb.ll2;
        dat.v.lb.type = 'll';
        
    elseif opt.optimise.v.v % dat.f.observed
    
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

        % Trace of inv(post)*prior
        % ------------------------
%         % Approximation where all off-diagonal elements of L are zero
%         dat.v.lb.tr = spm_diffeo('trapprox', h, double([opt.tpl.vs (model.mixreg.w(1)*model.v.l + model.mixreg.w(2))*opt.pg.prm]));
%         dat.v.lb.tr = dat.v.lb.tr(1);
        dat.v.lb.tr = trapprox((model.mixreg.w(1)*model.v.l + model.mixreg.w(2))*opt.pg.prm, h, 'vs', opt.tpl.vs);
        dat.v.lb.tr = dat.v.lb.tr / (model.mixreg.w(1)*model.v.l + model.mixreg.w(2));

        % LogDet of posterior precision
        % -----------------------------
%         % Approximation where all off-diagonal elements of L are zero
%         h(:,:,:,1) = h(:,:,:,1) * (model.mixreg.w(1)*model.v.l + model.mixreg.w(2)) * opt.pg.ker(1);
%         h(:,:,:,2) = h(:,:,:,2) * (model.mixreg.w(1)*model.v.l + model.mixreg.w(2)) * opt.pg.ker(2);
%         h(:,:,:,3) = h(:,:,:,3) * (model.mixreg.w(1)*model.v.l + model.mixreg.w(2)) * opt.pg.ker(3);
%         if opt.model.dim == 2
%             h(:,:,:,3) = 1;
%         end
%         h = spm_matcomp('Pointwise3', h, 'd');
%         h(h <= 0) = nan;
%         dat.v.lb.ld = sum(log(h(:)), 'omitnan');
        dat.v.lb.ld = ldapprox((model.mixreg.w(1)*model.v.l + model.mixreg.w(2))*opt.pg.prm, h, 'vs', opt.tpl.vs);
        clear h

        % Lower bound
        % -----------
        % KL divergence between multivariate normal distributions
        % posterior: mean = v , precision = H + (a1*l+a2)*L
        % prior:     a1 * [ mean = Wz , precision = l*L ]
        %            a2 * [ mean = 0  , precision = L   ]
        % (the prior is a mixture of two Gaussian distributions)
        K = prod(opt.tpl.lat)*3;
        dat.v.lb.uncty = trace(dat.z.S * model.pg.ww);
        if opt.optimise.v.l
            loglam = spm_prob('Gamma', 'Elog', model.v.l, model.v.n, K);
        else
            loglam = log(model.v.l);
        end
        dat.v.lb.ll1 = -0.5*( K*log(2*pi) ...
                              - K * loglam ...
                              - opt.pg.ld ...
                              + model.v.l * dat.v.lb.reg ...
                              + model.v.l * dat.v.lb.uncty ...
                              + model.v.l * dat.v.lb.tr );
        dat.v.lb.ll2 = -0.5*( K * log(2*pi) ...
                              - opt.pg.ld ...
                              + dat.v.lb.regv ...
                              + dat.v.lb.tr );
        dat.v.lb.val = -0.5*( - K ...
                              - model.mixreg.w(1) * K * loglam ...
                              - opt.pg.ld  ...
                              + dat.v.lb.ld ...
                              + model.mixreg.w(1) * model.v.l * dat.v.lb.reg ...
                              + model.mixreg.w(2) * dat.v.lb.regv ...
                              + model.mixreg.w(1) * model.v.l * dat.v.lb.tr ...
                              + model.mixreg.w(2) * dat.v.lb.tr ...
                              + model.mixreg.w(1) * model.v.l * dat.v.lb.uncty );
        dat.v.lb.type = 'kl';
    else
        dat.v.lb.val   = 0;
        dat.v.lb.tr    = 0;
        dat.v.lb.reg   = 0;
        dat.v.lb.uncty = 0;
        dat.v.lb.ll1   = 0;
        dat.v.lb.ll2   = 0;
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
        if opt.optimise.q.A
            logdetA = spm_prob('Wishart', 'ELogDet', model.q.A, model.q.n, 'normal');
        else
            logdetA = spm_matcomp('LogDet', model.q.A);
        end
        rind = opt.q.rind;
        qq   = dat.q.qq(rind,rind);
        Sq   = dat.q.S(rind,rind);
        dat.q.lb.type = 'kl';
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
    if isfield(model.lb, 'q'),    model.lb.q.val  = 0; end
    if isfield(model.lb, 'z'),    model.lb.z.val  = 0; end
    if isfield(model.lb, 'v1'),   model.lb.v1.val = 0; end
    if isfield(model.lb, 'v2'),   model.lb.v2.val = 0; end
    model.v.tr      = 0;
    model.v.uncty   = 0;
    model.v.reg     = 0;
    
    % --- Batch processing
    N = numel(dat);
    if opt.ui.verbose, before = plotBatchBegin('Init Lap'); end
    for i=1:ceil(N/batch)
        n1 = (i-1)*batch + 1;
        ne = min(N, i*batch);

        if opt.ui.verbose, before = plotBatch(i, batch, N, 50, before); end
    
        [opt.dist, dat(n1:ne)] = distribute(opt.dist, 'pgva_batch', 'oneInitLowerBound', 'inplace', dat(n1:ne), model, opt);
        
        for n=n1:ne
            
            % Add individual contributions
            % ----------------------------
            if isfield(model.lb, 'q'),    model.lb.q.val = model.lb.q.val + dat(n).q.lb.val;   end
            if isfield(model.lb, 'z'),    model.lb.z.val = model.lb.z.val + dat(n).z.lb.val;   end
            if isfield(model.v, 'uncty'), model.v.uncty  = model.v.uncty  + dat(n).v.lb.uncty; end
            if isfield(model.v, 'reg'),   model.v.reg    = model.v.reg    + dat(n).v.lb.reg;   end
            if dat(n).v.observed
                if isfield(model.lb, 'v2'),model.lb.v2.val = model.lb.v2.val + dat(n).v.lb.val; end
            else
                if isfield(model.v, 'tr'),  model.v.tr      = model.v.tr      + dat(n).v.lb.tr;  end
                if isfield(model.lb, 'v1'), model.lb.v1.val = model.lb.v1.val + dat(n).v.lb.val; end
            end
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
        end
    end
    if opt.ui.verbose, plotBatchEnd; end

end

%% ------------------------------------------------------------------------
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
    
    % -------------------------
    % CASE 1: observed velocity
    % -------------------------
    if dat.v.observed
        
        % Lower bound
        % -----------
        % Log-likelihood of a multivariate normal distribution
        % mean = W*z, precision = l*L
        if opt.pg.provided
            wz = reconstructVelocity('latent', dat.z.z, 'subspace', model.pg.w, ...
                                     'loop', loop, 'par', par);
            v = numeric(dat.v.v); 
            if model.mixreg.w(2)
                m = numeric(dat.v.m);
                dat.v.lb.regv = v(:)' * m(:);
            else
                dat.v.lb.regv = 0;
            end
            r = v - wz;
            clear v wz
            m = spm_diffeo('vel2mom', r, double([opt.tpl.vs opt.pg.prm]));
            dat.v.lb.reg = r(:)' * m(:);
            clear r m
        else
            dat.v.lb.reg  = single(dat.v.v(:))' * single(dat.v.m(:));
            dat.v.lb.regv = 0;
        end
    
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
        dat.v.ok  = 1; % for GN failure tracking
        dat.v.ok2 = 0; % for GN failure tracking
        dat.v.v = prepareOnDisk(dat.v.v, [opt.tpl.lat 3]);
        dat.v.m = prepareOnDisk(dat.v.m, [opt.tpl.lat 3]);
        if strcmpi(mode, 'zero')
            dat.v.v(:)    = 0;
            dat.v.m(:)    = 0;
            dat.v.lb.reg  = 0;
            dat.v.lb.regv = 0;
        else
            v = create([opt.tpl.lat 3], 'single');
            dat.v.v(:,:,:,:) = v;
            m = spm_diffeo('vel2mom', v, double([opt.tpl.vs opt.pg.prm]));
            dat.v.m(:,:,:,:) = m;
            if model.mixreg.w(2)
                dat.v.lb.regv = v(:)' * m(:);
            else
                dat.v.lb.regv = 0;
            end
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
    
    end

end

function [dat, model] = batchInitVelocity(mode, dat, model, opt)

    % --- Detect parallelisation scheme
    if strcmpi(opt.split.loop, 'subject') && opt.split.par > 0
        batch = opt.split.batch;
    else
        batch = 1;
    end
    
    % --- Batch processing
    N = numel(dat);
    if opt.ui.verbose, before = plotBatchBegin('Init Vel'); end
    for i=1:ceil(N/batch)
        n1 = (i-1)*batch + 1;
        ne = min(N, i*batch);

        if opt.ui.verbose, before = plotBatch(i, batch, N, 50, before); end
    
        [opt.dist, dat(n1:ne)] = distribute(opt.dist, 'pgva_batch', 'OneInitVelocity', 'inplace', dat(n1:ne), model, opt, mode);
        
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

%% ========================================================================
%    E-step
% =========================================================================

%% ------------------------------------------------------------------------
% Fit Affine
% -------------------------------------------------------------------------

function dat = oneFitAffine(dat, model, opt)

    % If v is observed
    % ----------------
    if dat.v.observed
        return
    end
    
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
    if opt.iter.pena && dat.q.ok < 0
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
            if opt.optimise.q.A
                logdetA = spm_prob('Wishart', 'ELogDet', model.q.A, model.q.n, 'normal');
            else
                logdetA = spm_matcomp('LogDet', model.q.A);
            end
            rind = opt.q.rind;
            qq   = dat.q.qq(rind,rind);
            Sq   = dat.q.S(rind,rind);
            dat.q.lb.type = 'kl';
            dat.q.lb.val  = -0.5*( trace((Sq + qq) * model.q.A) ...
                                   - logdetA ...
                                   + spm_matcomp('LogDet', Sq) ...
                                   - opt.q.Mr );
        end
                              
    end % < penalise previous failure
    
    % Cleaning (just in case)
    % --------
    dat.v.iphi = rmarray(dat.v.iphi);
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
        [opt.dist, dat(n1:ne)] = distribute(opt.dist, 'pgva_batch', 'oneFitAffine', 'inplace', dat(n1:ne), model, opt);
        
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
    
    % Update covariance
    % -----------------
    dat.z.S = model.mixreg.w(1) * model.v.l * model.pg.ww + model.z.A;
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
    dat.z.z  = model.mixreg.w(1) * model.v.l * dat.z.S * wm;
    dat.z.zz = dat.z.z * dat.z.z';
    
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
        [opt.dist, dat(n1:ne)] = distribute(opt.dist, 'pgva_batch', 'oneFitLatent', 'inplace', dat(n1:ne), model, opt);
        
        for n=n1:ne
            
            % Add individual contributions
            % ----------------------------
            model.z.z      = model.z.z      + dat(n).z.z;
            model.z.zz     = model.z.zz     + dat(n).z.zz;
            model.z.S      = model.z.S      + dat(n).z.S;
            model.z.Z(:,n) = dat(n).z.z;
            
        end
        
    end
    if opt.ui.verbose, plotBatchEnd; end
end


%% ------------------------------------------------------------------------
% Fit Velocity field
% -------------------------------------------------------------------------

function dat = oneFitVelocity(dat, model, opt)

    % If v is observed
    % ----------------
    if dat.v.observed
        return
    end
    
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
    if opt.iter.pena && dat.v.ok < 0
        dat.v.ok = dat.v.ok + 1;
    else
        
        % Gauss-Newton iterations
        % -----------------------
        cumok = false;
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

            if checkarray(model.pg.w) && all(dat.z.z ~= 0)
                wz = reconstructVelocity('latent', dat.z.z, 'subspace', model.pg.w, ...
                                         'loop', loop, 'par', par);
                g = g - spm_diffeo('vel2mom', single(wz), double([opt.tpl.vs (model.mixreg.w(1) * model.v.l * opt.pg.prm)]));
                v = numeric(dat.v.v);
                r = v - wz;
                clear wz
            else
                r = dat.v.v;
                v = dat.v.v;
            end
            g = g + ghPriorVel(v, opt.tpl.vs, (model.mixreg.w(1)*model.v.l + model.mixreg.w(2)) * opt.pg.prm, opt.pg.bnd);

            % Compute search direction
            % ------------------------
            dv = -spm_diffeo('fmg', single(h), single(g), ...
                double([opt.tpl.vs (model.mixreg.w(1)*model.v.l + model.mixreg.w(2)) * opt.pg.prm 2 2]));
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
                'v0', v, 'lam', model.mixreg.w(1) * model.v.l, 'geod', model.mixreg.w(2), ...
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
                m = spm_diffeo('vel2mom', result.v, double([opt.tpl.vs opt.pg.prm]));
                if model.mixreg.w(2)
                    dat.v.lb.regv = result.v(:)' * m(:);
                end
                dat.v.m       = copyarray(m, dat.v.m);
                clear m
                dat.v.r       = copyarray(result.r,    dat.v.r);
                dat.v.ipsi    = copyarray(result.ipsi, dat.v.ipsi);
                dat.f.pf      = copyarray(result.pf,   dat.f.pf);
                dat.f.c       = copyarray(result.c,    dat.f.c);
                dat.f.bb      = result.bb;
                if strcmpi(opt.match, 'pull')
                    dat.tpl.wmu   = copyarray(result.wmu, dat.tpl.wmu);
                    rmarray(result.wa);
                end
                r = result.r;
                ipsi = dat.v.ipsi;
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
    
        % -----------
        % Lower bound
        % -----------
        % KL divergence between multivariate normal distributions
        % posterior: mean = v,    precision = H + l*L
        % prior:     mean = W*z , precision = l*L
        
        % Regularisation part
        % -------------------
        m = spm_diffeo('vel2mom', single(numeric(r)), double([opt.tpl.vs opt.pg.prm]));
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
%         dat.v.lb.tr = spm_diffeo('trapprox', h, double([opt.tpl.vs (model.mixreg.w(1)*model.v.l + model.mixreg.w(2)) * opt.pg.prm]));
%         dat.v.lb.tr = dat.v.lb.tr(1);
        dat.v.lb.tr = trapprox((model.mixreg.w(1)*model.v.l + model.mixreg.w(2))*opt.pg.prm, h, 'vs', opt.tpl.vs);
        dat.v.lb.tr = dat.v.lb.tr / (model.mixreg.w(1)*model.v.l + model.mixreg.w(2));
        
        % LogDet(P)
        % ---------
%         % Approximation where all off-diagonal elements of L are zero
%         h(:,:,:,1) = h(:,:,:,1) * (model.mixreg.w(1)*model.v.l + model.mixreg.w(2)) * opt.pg.ker(1);
%         h(:,:,:,2) = h(:,:,:,2) * (model.mixreg.w(1)*model.v.l + model.mixreg.w(2)) * opt.pg.ker(2);
%         h(:,:,:,3) = h(:,:,:,3) * (model.mixreg.w(1)*model.v.l + model.mixreg.w(2)) * opt.pg.ker(3);
%         if opt.model.dim == 2
%             h(:,:,:,3) = 1;
%         end
%         h = spm_matcomp('Pointwise3', h, 'd');
%         h(h <= 0) = nan;
%         dat.v.lb.ld = sum(log(h(:)), 'omitnan');
        dat.v.lb.ld = ldapprox((model.mixreg.w(1)*model.v.l + model.mixreg.w(2))*opt.pg.prm, h, 'vs', opt.tpl.vs);
        clear h
        
        % KL divergence
        % -------------
        K = prod(opt.tpl.lat) * 3;
        if opt.optimise.v.l
            loglam = spm_prob('Gamma', 'Elog', model.v.l, model.v.n, K);
        else
            loglam = log(model.v.l);
        end
        dat.v.lb.ll1 = -0.5*( K*log(2*pi) ...
                              - K *loglam ...
                              - opt.pg.ld ...
                              + model.v.l * dat.v.lb.reg ...
                              + model.v.l * dat.v.lb.uncty ...
                              + model.v.l * dat.v.lb.tr );
        dat.v.lb.ll2 = -0.5*( K * log(2*pi) ...
                              - opt.pg.ld ...
                              + dat.v.lb.regv ...
                              + dat.v.lb.tr );
        dat.v.lb.val = -0.5*( - K ...
                              - model.mixreg.w(1) * K * loglam ...
                              - opt.pg.ld  ...
                              + dat.v.lb.ld ...
                              + model.mixreg.w(1) * model.v.l * dat.v.lb.reg ...
                              + model.mixreg.w(2) * dat.v.lb.regv ...
                              + model.mixreg.w(1) * model.v.l * dat.v.lb.tr ...
                              + model.mixreg.w(2) * dat.v.lb.tr ...
                              + model.mixreg.w(1) * model.v.l * dat.v.lb.uncty );
    
    end % < penalise previous failure
    
    % Cleaning
    % --------
    % Just in case
    dat.v.iphi = rmarray(dat.v.iphi);
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
        [opt.dist, dat(n1:ne)] = distribute(opt.dist, 'pgva_batch', 'oneFitVelocity', 'inplace', dat(n1:ne), model, opt);
        
        
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
            if dat.f.observed
                if strcmpi(opt.match, 'pull')
                    dat.f.lb.val = llMatching(noisemodel, dat.tpl.wmu, dat.f.f, ...
                        'par', par, 'loop', loop, 'debug', opt.ui.debug);
                else
                    dat.f.lb.val = llMatching(noisemodel, model.tpl.mu, dat.f.pf, dat.f.c, ...
                        'bb', dat.f.bb, 'par', par, 'loop', loop, 'debug', opt.ui.debug);
                end
            end
            
        case 'template'
            if strcmpi(opt.match, 'pull')
                if opt.tpl.cat
                    dat.tpl.wa = pullTemplate(dat.v.ipsi, model.tpl.a, ...
                        'par', par, 'output', dat.tpl.wa, 'debug', opt.ui.debug);
                    [dat.tpl.wmu,ssmu] = reconstructProbaTemplate(dat.tpl.wa, ...
                        'scale', dat.tpl.scale, ...
                        'loop', loop, 'par', par, 'output', dat.tpl.wmu, ...
                        'debug', opt.ui.debug);
                    if opt.optimise.tpl.scale
                        % Individual probability scaling
                        % ------------------------------
                        for i=1:opt.iter.scl
                            dimf = [size(dat.f.f) 1];
                            dat.tpl.scale = (dat.f.ss1(:) * dat.f.ss0) ./ (ssmu(:) * prod(dimf(1:3)));
                            [dat.tpl.wmu, ssmu] = reconstructProbaTemplate(dat.tpl.wa, ...
                                'scale', dat.tpl.scale, ...
                                'loop', loop, 'par', par, 'output', dat.tpl.wmu, ...
                                'debug', opt.ui.debug);
                        end
                    end
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
                if opt.optimise.q.A
                    logdetA = spm_prob('Wishart', 'ELogDet', model.q.A, model.q.n, 'normal');
                else
                    logdetA = spm_matcomp('LogDet', model.q.A);
                end
                rind = opt.q.rind;
                qq   = dat.q.qq(rind,rind);
                Sq   = dat.q.S(rind,rind);
                dat.q.lb.type = 'kl';
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
            wz = reconstructVelocity('latent', dat.z.z, 'subspace', model.pg.w, ...
                                     'loop', loop, 'par', par);
            r = numeric(dat.v.v) - wz;
            clear wz
            m = spm_diffeo('vel2mom', r, double([opt.tpl.vs opt.pg.prm]));
            dat.v.lb.reg = r(:)' * m(:);
            clear r m
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
            K = prod(opt.tpl.lat) * 3;
            if opt.optimise.v.l
                loglam = spm_prob('Gamma', 'Elog', model.v.l, model.v.n, K);
            else
                loglam = log(model.v.l);
            end
            if dat.v.observed
                dat.v.lb.ll1 = -0.5*( K*log(2*pi) ...
                                      - K * loglam ...
                                      - opt.pg.ld ...
                                      + model.v.l * dat.v.lb.reg ...
                                      + model.v.l * dat.v.lb.uncty );
                dat.v.lb.ll2 = -0.5*( K * log(2*pi) ...
                                      - opt.pg.ld ...
                                      + dat.v.lb.regv );
                dat.v.lb.val = model.mixreg.w(1) * dat.v.lb.ll1 + ...
                               model.mixreg.w(2) * dat.v.lb.ll2;
            else
                dat.v.lb.ll1 = -0.5*( K*log(2*pi) ...
                                      - K * loglam ...
                                      - opt.pg.ld ...
                                      + model.v.l * dat.v.lb.reg ...
                                      + model.v.l * dat.v.lb.uncty ...
                                      + model.v.l * dat.v.lb.tr );
                dat.v.lb.ll2 = -0.5*( K * log(2*pi) ...
                                      - opt.pg.ld ...
                                      + dat.v.lb.regv ...
                                      + dat.v.lb.tr );
                dat.v.lb.val = -0.5*( - K ...
                                      - model.mixreg.w(1) * K * loglam ...
                                      - opt.pg.ld  ...
                                      + dat.v.lb.ld ...
                                      + model.mixreg.w(1) * model.v.l * dat.v.lb.reg ...
                                      + model.mixreg.w(2) * dat.v.lb.regv ...
                                      + model.mixreg.w(1) * model.v.l * dat.v.lb.tr ...
                                      + model.mixreg.w(2) * dat.v.lb.tr ...
                                      + model.mixreg.w(1) * model.v.l * dat.v.lb.uncty );
            end
            
        case 'latent'
            wz = reconstructVelocity('latent', dat.z.z, 'subspace', model.pg.w, ...
                                     'loop', loop, 'par', par);
            r = numeric(dat.v.v) - wz;
            clear wz
            m = spm_diffeo('vel2mom', r, double([opt.tpl.vs opt.pg.prm]));
            dat.v.lb.reg = r(:)' * m(:);
            clear r m
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
            K = prod(opt.tpl.lat) * 3;
            if opt.optimise.v.l
                loglam = spm_prob('Gamma', 'Elog', model.v.l, model.v.n, K);
            else
                loglam = log(model.v.l);
            end
            if dat.v.observed
                dat.v.lb.ll1 = -0.5*( K*log(2*pi) ...
                                      - K * loglam ...
                                      - opt.pg.ld ...
                                      + model.v.l * dat.v.lb.reg ...
                                      + model.v.l * dat.v.lb.uncty );
                dat.v.lb.ll2 = -0.5*( K * log(2*pi) ...
                                      - opt.pg.ld ...
                                      + dat.v.lb.regv );
                dat.v.lb.val = model.mixreg.w(1) * dat.v.lb.ll1 + ...
                               model.mixreg.w(2) * dat.v.lb.ll2;
            else
                dat.v.lb.ll1 = -0.5*( K*log(2*pi) ...
                                      - K * loglam ...
                                      - opt.pg.ld ...
                                      + model.v.l * dat.v.lb.reg ...
                                      + model.v.l * dat.v.lb.uncty ...
                                      + model.v.l * dat.v.lb.tr );
                dat.v.lb.ll2 = -0.5*( K * log(2*pi) ...
                                      - opt.pg.ld ...
                                      + dat.v.lb.regv ...
                                      + dat.v.lb.tr );
                dat.v.lb.val = -0.5*( - K ...
                                      - model.mixreg.w(1) * K * loglam ...
                                      - opt.pg.ld  ...
                                      + dat.v.lb.ld ...
                                      + model.mixreg.w(1) * model.v.l * dat.v.lb.reg ...
                                      + model.mixreg.w(2) * dat.v.lb.regv ...
                                      + model.mixreg.w(1) * model.v.l * dat.v.lb.tr ...
                                      + model.mixreg.w(2) * dat.v.lb.tr ...
                                      + model.mixreg.w(1) * model.v.l * dat.v.lb.uncty );
            end
            if opt.optimise.z.A
                logdetA = spm_prob('Wishart', 'ELogDet', model.z.A, model.z.n, 'normal');
            else
                logdetA = spm_matcomp('LogDet', model.z.A);
            end
            dat.z.lb.val  = -0.5*( trace((dat.z.S + dat.z.zz) * model.z.A) ...
                                   - logdetA ...
                                   - spm_matcomp('LogDet', dat.z.S) ...
                                   - opt.pg.K );
            
        case {'lambda', 'mixture'}
            K = prod(opt.tpl.lat) * 3;
            if opt.optimise.v.l
                loglam = spm_prob('Gamma', 'Elog', model.v.l, model.v.n, K);
            else
                loglam = log(model.v.l);
            end
            if dat.v.observed
                dat.v.lb.ll1 = -0.5*( K*log(2*pi) ...
                                      - K * loglam ...
                                      - opt.pg.ld ...
                                      + model.v.l * dat.v.lb.reg ...
                                      + model.v.l * dat.v.lb.uncty );
                dat.v.lb.ll2 = -0.5*( K * log(2*pi) ...
                                      - opt.pg.ld ...
                                      + dat.v.lb.regv );
                dat.v.lb.val = model.mixreg.w(1) * dat.v.lb.ll1 + ...
                               model.mixreg.w(2) * dat.v.lb.ll2;
            else
                dat.v.lb.ll1 = -0.5*( K*log(2*pi) ...
                                      - K * loglam ...
                                      - opt.pg.ld ...
                                      + model.v.l * dat.v.lb.reg ...
                                      + model.v.l * dat.v.lb.uncty ...
                                      + model.v.l * dat.v.lb.tr );
                dat.v.lb.ll2 = -0.5*( K * log(2*pi) ...
                                      - opt.pg.ld ...
                                      + dat.v.lb.regv ...
                                      + dat.v.lb.tr );
                dat.v.lb.val = -0.5*( - K ...
                                      - model.mixreg.w(1) * K * loglam ...
                                      - opt.pg.ld  ...
                                      + dat.v.lb.ld ...
                                      + model.mixreg.w(1) * model.v.l * dat.v.lb.reg ...
                                      + model.mixreg.w(2) * dat.v.lb.regv ...
                                      + model.mixreg.w(1) * model.v.l * dat.v.lb.tr ...
                                      + model.mixreg.w(2) * dat.v.lb.tr ...
                                      + model.mixreg.w(1) * model.v.l * dat.v.lb.uncty );
            end
            
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
    if isfield(model.lb, 'm'),    model.lb.m.val  = 0; end
    if isfield(model.lb, 'z'),    model.lb.z.val  = 0; end
    if isfield(model.lb, 'q'),    model.lb.q.val  = 0; end
    if isfield(model.lb, 'v1'),   model.lb.v1.val = 0; end
    if isfield(model.lb, 'v2'),   model.lb.v2.val = 0; end
    if isfield(model.v, 'tr'),    model.v.tr      = 0; end
    if isfield(model.v, 'reg'),   model.v.reg     = 0; end
    if isfield(model.v, 'uncty'), model.v.uncty   = 0; end
    
    N = numel(dat);
    
    % --- Batch processing
    if opt.ui.verbose, before = plotBatchBegin(['LB ' which(1:min(numel(which),3))]); end
    for i=1:ceil(N/batch)
        n1 = (i-1)*batch + 1;
        ne = min(N, i*batch);

        if opt.ui.verbose, before = plotBatch(i, batch, N, 50, before); end
    
        % Compute subjects grad/hess w.r.t. initial velocity
        % --------------------------------------------------
        [opt.dist, dat(n1:ne)] = distribute(opt.dist, 'pgva_batch', 'oneLB', 'inplace', dat(n1:ne), model, opt, which);
        
        for n=n1:ne
            
            % Add individual contributions
            % ----------------------------
            if isfield(model.lb, 'z'),    model.lb.z.val  = model.lb.z.val + dat(n).z.lb.val;   end
            if isfield(model.v, 'reg'),   model.v.reg     = model.v.reg    + dat(n).v.lb.reg;   end
            if isfield(model.v, 'uncty'), model.v.uncty   = model.v.uncty  + dat(n).v.lb.uncty; end
            if dat(n).v.observed
                if isfield(model.lb, 'v2'), model.lb.v2.val = model.lb.v2.val + dat(n).v.lb.val; end
            else % dat(n).f.observed
                if isfield(model.lb, 'm'),  model.lb.m.val  = model.lb.m.val  + dat(n).f.lb.val; end
                if isfield(model.lb, 'v1'), model.lb.v1.val = model.lb.v1.val + dat(n).v.lb.val; end
                if isfield(model.v, 'tr'),  model.v.tr      = model.v.tr      + dat(n).v.lb.tr;  end
                if isfield(model.lb, 'q'),  model.lb.q.val  = model.lb.q.val + dat(n).q.lb.val;  end
            end
            
        end
        
    end
    if opt.ui.verbose, plotBatchEnd; end

    % Model specific parts
    % --------------------
    switch lower(which)
        case 'orthogonalise'
            model.lb.w.val = llPriorSubspace(model.pg.w, model.pg.n * model.pg.ww, opt.pg.ld + prod(opt.tpl.lat)*3*log(model.pg.n));
            if opt.optimise.z.A && opt.z.n0
                model.lb.Az.val = -spm_prob('Wishart', 'kl', ...
                                            model.z.A,   model.z.n, ...
                                            opt.z.A0,    opt.z.n0, ...
                                            'normal');
            end
        case 'subspace'
            model.lb.w.val = llPriorSubspace(model.pg.w, model.pg.n * model.pg.ww, opt.pg.ld + prod(opt.tpl.lat)*3*log(model.pg.n));
        case 'precisionz'
            if opt.optimise.z.A && opt.z.n0
                model.lb.Az.val = -spm_prob('Wishart', 'kl', ...
                                            model.z.A,   model.z.n, ...
                                            opt.z.A0,    opt.z.n0, ...
                                            'normal');
            end
        case 'precisionq'
            if opt.optimise.q.A && opt.q.n0 && opt.q.Mr
                model.lb.Aq.val = -spm_prob('Wishart', 'kl', ...
                                            model.q.A,   model.q.n, ...
                                            opt.q.A0,    opt.q.n0, ...
                                            'normal');
            end
        case 'lambda'
            if opt.optimise.v.l && opt.v.n0
                model.lb.l.val  = -spm_prob('Gamma', 'kl', ...
                                            model.v.l, model.v.n, ...
                                            opt.v.l0,  opt.v.n0, ...
                                            prod(opt.tpl.lat)*3, 'normal');
            end
    end
end
