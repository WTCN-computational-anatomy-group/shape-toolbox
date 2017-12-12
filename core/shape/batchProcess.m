function varargout = batchProcess(id, varargin)
% _________________________________________________________________________
%
% Collection of functions that require updating some subject-specific
% variables. Depending on the needs, as much as possible is done in order
% to minimise i/o, memory and computational costs.
%
% FORMAT [dat, model] = batchProcess('InitAffine', mode, dat, model, opt)
%   Initialise affine coordinates (zero or randomly)
%   Update>> dat.q, dat.qq, dat.Sq, model.q, model.qq, model.Sq
%
% FORMAT [dat, model] = batchProcess('InitResidual', mode, dat, model, opt)
%   Initialise residual fields (zero or randomly)
%   Update>> dat.r, model.llr
%
% FORMAT [dat, model] = batchProcess('InitLatent', mode, dat, model, opt)
%   Initialise latent coordinates and insure they are zero-centered (but
%   orthogonalisation is not done here).
%   Update>> dat.z, dat.zz, dat.Sz, model.z, model.zz, model.Sz
%
% FORMAT dat = batchProcess('RotateLatent', dat, opt, R)
%   Apply a rotation matrix to all latent coordinates
%   Update>> dat.z, dat.zz, dat.Sz
%
% FORMAT dat = batchProcess('CentreLatent', dat, opt, (mean))
%   Subtract the mean to all latent coordinates
%   Update>> dat.z, dat.zz
%
% FORMAT [dat, model] = batchProcess('FitAffine', dat, model, opt)
%   Gauss-Newton update of the affine coordinates
%   Update>> dat.okq, dat.q, dat.qq, dat.Sq, dat.A, dat.llm, dat.pf, dat.c
%            model.q, model.qq, model.Sq, model.llm
%
% FORMAT [dat, model] = batchProcess('FitLatent', dat, model, opt)
%   Gauss-Newton update of the latent coordinates
%   Update>> dat.okz, dat.z, dat.zz, dat.Sz, dat.llm, dat.pf, dat.c
%            model.z, model.zz, model.Sz, model.llm
%
% FORMAT [dat, model] = batchProcess('FitResidual', dat, model, opt)
%   Gauss-Newton update of the residual fields
%   Update>> dat.okr, dat.r, dat.llr, dat.llm, dat.pf, dat.c
%            model.llm, model.llr
%
% FORMAT [dat, model] = batchProcess('GradHessSubspace', dat, model, opt)
%   Gradient and Hessian of the matching part w.r.t. subspace
%   Update>> model.gw, model.hw
%
% FORMAT dat = batchProcess('Update', dat, model, opt, {...}, 'clean', {...})
%   Generic tool to update specified subject variables. Names of variables
%   to update are provided in the first list, and names of variables to
%   discard in the end (to save memory/disk) can be provided in the second
%   list.
% _________________________________________________________________________
%
% The following subfunctions should not be used usually. They are only
% needed when distributing jobs on a cluster.
% 
% FORMAT dat = batchProcess('OneStepFitAffine', dat, model, opt)
% FORMAT dat = batchProcess('OneStepFitLatent', dat, model, opt)
% FORMAT dat = batchProcess('OneStepFitResidual', dat, model, opt)
% FORMAT dat = batchProcess('OneStepGradHessVelocity', dat, model, opt)
% FORMAT dat = batchProcess('OneUpdate', dat, model, opt, todo, toclean)
% _________________________________________________________________________


    switch lower(id)
        % -------
        %  BATCH
        % -------
        case 'initaffine'
            [varargout{1:nargout}] = batchInitAffine(varargin{:});
        case 'initlatent'
            [varargout{1:nargout}] = batchInitLatent(varargin{:});
        case 'initresidual'
            [varargout{1:nargout}] = batchInitResidual(varargin{:});
        case 'initzero'
            [varargout{1:nargout}] = batchInitZero(varargin{:});
        case 'rotatelatent'
            [varargout{1:nargout}] = batchRotateLatent(varargin{:});
        case 'centrelatent'
            [varargout{1:nargout}] = batchCentreLatent(varargin{:});
        case 'fitaffine'
            [varargout{1:nargout}] = batchFitAffine(varargin{:});
        case 'fitlatent'
            [varargout{1:nargout}] = batchFitLatent(varargin{:});
        case 'fitresidual'
            [varargout{1:nargout}] = batchFitResidual(varargin{:});
        case 'gradhesssubspace'
            [varargout{1:nargout}] = batchGradHessSubspace(varargin{:});
        case 'update'
            [varargout{1:nargout}] = batchUpdate(varargin{:});
        % ------
        %  STEP
        % ------
        case 'onestepinitresidual'
            [varargout{1:nargout}] = oneStepInitResidual(varargin{:});
        case 'onestepinitzero'
            [varargout{1:nargout}] = oneStepInitZero(varargin{:});
        case 'onestepfitaffine'
            [varargout{1:nargout}] = oneStepFitAffine(varargin{:});
        case 'onestepfitlatent'
            [varargout{1:nargout}] = oneStepFitLatent(varargin{:});
        case 'onestepfitresidual'
            [varargout{1:nargout}] = oneStepFitResidual(varargin{:});
        case 'onestepgradhessvelocity'
            [varargout{1:nargout}] = oneStepGradHessVelocity(varargin{:});
        case 'oneupdate'
            [varargout{1:nargout}] = oneUpdate(varargin{:});
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
    fprintf(' | %fs\n', toc);
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
            error('Unknwon mode %s', mode)
    end
    
    % Init model
    % ----------
    M        = size(opt.affine_basis, 3);
    model.q  = zeros(M, 1);
    model.qq = zeros(M);
    model.Sq = zeros(M);
    
    % Init subjects
    % -------------
    if opt.verbose, before = plotBatchBegin('Init Q'); end;
    for n=1:opt.N
        if opt.verbose, before = plotBatch(n, 1, opt.N, 50, before); end;
        q         = create([M, 1]);
        dat(n).q  = saveOnDisk(dat(n).q,  q);
        dat(n).qq = saveOnDisk(dat(n).qq, q*q');
        dat(n).Sq = saveOnDisk(dat(n).Sq, zeros(M));
        model.q   = model.q  + q;
        model.qq  = model.qq + q*q';
        dat(n).A  = exponentiateAffine(q, opt.affine_basis);
    end
    if opt.verbose, plotBatchEnd; end;

end

% -------------------------------------------------------------------------
%    Init Residual
% -------------------------------------------------------------------------

function dat = oneStepInitResidual(dat, ~, opt, mode)

    % --- Don't parallelise here, there is no need
    
    switch lower(mode)
        case 'zero'
            create = @zeros;
        case 'rand'
            create = @randn;
        otherwise
            error('Unknwon mode %s', mode)
    end
    
    if strcmpi(mode, 'zero')
        dat.r = prepareOnDisk(dat.r, [opt.lat 3]);
        dat.r(:) = 0;
    else
        r         = create([opt.lat 3]);
        dat.r  = saveOnDisk(dat.r, r);
    end
%     if isfield(opt, 'logdet')
%         dat.llr = llPriorVelocity(r, ...
%             'vs', opt.vs,  'logdet', opt.logdet, ...
%             'prm', opt.prm, 'debug', opt.debug);
%     else
%         dat.llr = llPriorVelocity(r, ...
%             'vs', opt.vs, 'prm', opt.prm, 'debug', opt.debug);
%     end

end

function [dat, model] = batchInitResidual(mode, dat, model, opt)
    % --- Detect parallelisation scheme
    if strcmpi(opt.loop, 'subject') && opt.par > 0
        batch = opt.batch;
    else
        batch = 1;
    end

%     % Init log-likelihood
%     % -------------------
%     model.llr = 0;
    
    % --- Batch processing
    if opt.verbose, before = plotBatchBegin('Init R'); end;
    for i=1:ceil(opt.N/batch)
        n1 = (i-1)*batch + 1;
        ne = min(opt.N, i*batch);

        if opt.verbose, before = plotBatch(i, batch, opt.N, 50, before); end;
    
        % Compute subjects grad/hess w.r.t. initial velocity
        % --------------------------------------------------
        dat(n1:ne) = distribute('OneStepInitResidual', dat(n1:ne), model, opt, mode);
        
%         for n=n1:ne
%             
%             % Add individual contributions
%             % ----------------------------
%             model.llr = model.llr + dat(n).llr;
%         end
        
    end
    if opt.verbose, plotBatchEnd; end;
end

% -------------------------------------------------------------------------
%    Init Zeros
% -------------------------------------------------------------------------
% Init        r    = 0
%             v    = 0
%             iphi = id
% Reconstruct ipsi
%             pf c
% Compute     llm
% Clean       iphi
%             ipsi

function dat = oneStepInitZero(dat, model, opt)
    
    if strcmpi(opt.loop, 'subject')
        loop = '';
        par  = 0;
    else
        loop = opt.loop;
        par  = opt.par;
    end
    
    dat.r    = prepareOnDisk(dat.r, [opt.lat 3]);
    dat.r(:) = 0;
    dat.v    = prepareOnDisk(dat.v, [opt.lat 3]);
    dat.v(:) = 0;
    iphi     = spm_warps('identity', opt.lat);
    latf = [size(dat.f) 1];
    latf = latf(1:3);
    ipsi = reconstructIPsi(dat.A, iphi, ...
        'lat', latf, 'Mf', dat.Mf, 'Mmu', model.Mmu, ...
        'debug', opt.debug);
    clear iphi
    [dat.pf, dat.c, dat.bb] = pushImage(ipsi, dat.f, opt.lat, ...
        'loop', loop, 'par', par, ...
        'output', {dat.pf, dat.c}, 'debug', opt.debug);

end

function [dat, model] = batchInitZero(dat, model, opt)
    % --- Detect parallelisation scheme
    if strcmpi(opt.loop, 'subject') && opt.par > 0
        batch = opt.batch;
    else
        batch = 1;
    end
    
    % --- Batch processing
    if opt.verbose, before = plotBatchBegin('Init Zero'); end;
    for i=1:ceil(opt.N/batch)
        n1 = (i-1)*batch + 1;
        ne = min(opt.N, i*batch);

        if opt.verbose, before = plotBatch(i, batch, opt.N, 50, before); end;
    
        % Compute subjects grad/hess w.r.t. initial velocity
        % --------------------------------------------------
        dat(n1:ne) = distribute('OneStepInitZero', dat(n1:ne), model, opt);
        
    end
    if opt.verbose, plotBatchEnd; end;
end

% -------------------------------------------------------------------------
%    Init Latent
% -------------------------------------------------------------------------

function [dat, model] = batchInitLatent(mode, dat, model, opt)

    % --- Don't parallelise here, there is no need
    
    switch lower(mode)
        case 'zero'
            create = @zeros;
        case 'rand'
            create = @randn;
        otherwise
            error('Unknwon mode %s', mode)
    end
    
    % Init model
    % ----------
    mz  = zeros(opt.K, 1);
    mzz = zeros(opt.K);
    
    % Init subjects
    % -------------
    if opt.verbose, before = plotBatchBegin('Init Z'); end;
    for n=1:opt.N
        if opt.verbose, before = plotBatch(n, 1, opt.N, 50, before); end;
        z        = create([opt.K, 1]);
        dat(n).z = saveOnDisk(dat(n).z, z);
        mz       = mz  + z;
        mzz      = mzz + z*z';
    end
    if opt.verbose, plotBatchEnd; end;
    
    % Center subjects
    % ---------------
    if opt.verbose, before = plotBatchBegin('Center Z'); end;
    for n=1:opt.N
        if opt.verbose, before = plotBatch(n, 1, opt.N, 50, before); end;
        z = numeric(dat(n).z) - mz/opt.N;
        
        dat(n).z  = saveOnDisk(dat(n).z, z);
        dat(n).zz = saveOnDisk(dat(n).z, z*z');
        dat(n).Sz = saveOnDisk(dat(n).Sz, zeros(opt.K));
    end
    if opt.verbose, plotBatchEnd; end;
    model.z  = saveOnDisk(model.z, zeros(opt.K, 1));
    model.zz = saveOnDisk(model.zz, mzz - mz*mz'/opt.N);
    model.Sz = saveOnDisk(model.Sz, zeros(opt.K));

end

% -------------------------------------------------------------------------
%    Center
% -------------------------------------------------------------------------

function dat = batchCentreLatent(dat, opt, mean)

    % --- Don't parallelise here, there is no need
    
    % Compute the mean if needed
    % --------------------------
    if nargin < 3
        if opt.verbose, before = plotBatchBegin('Accumul Z'); end;
        mean = zeros(opt.K, 1);
        for n=1:opt.N
            if opt.verbose, before = plotBatch(n, 1, opt.N, 50, before); end;
            mean = mean + numeric(dat(n).z);
        end
        if opt.verbose, plotBatchEnd; end;
        mean = mean / opt.N;
    end
        
    
    % Center subjects
    % ---------------
    if opt.verbose, before = plotBatchBegin('Centre Z'); end;
    for n=1:opt.N
        if opt.verbose, before = plotBatch(n, 1, opt.N, 50, before); end;
        z = numeric(dat(n).z) - mean;
        dat(n).z  = saveOnDisk(dat(n).z, z);
        dat(n).zz = saveOnDisk(dat(n).z, z*z');
    end
    if opt.verbose, plotBatchEnd; end;
end

% -------------------------------------------------------------------------
%    Rotate
% -------------------------------------------------------------------------

function dat = batchRotateLatent(dat, opt, R)

    % --- Don't parallelise here, there is no need
    
    % Rotate subjects
    % ---------------
    if opt.verbose, before = plotBatchBegin('Rotate Z'); end;
    for n=1:opt.N
        if opt.verbose, before = plotBatch(n, 1, opt.N, 50, before); end;
        z = R * numeric(dat(n).z);
        dat(n).z  = saveOnDisk(dat(n).z, z);
        dat(n).zz = saveOnDisk(dat(n).z, z*z');
        dat(n).Sz = saveOnDisk(dat(n).Sz, R * numeric(dat(n).Sz) * R');
    end
    if opt.verbose, plotBatchEnd; end;
end

% -------------------------------------------------------------------------
%    E-step
% -------------------------------------------------------------------------

% ------
% Affine
% ------

function dat = oneStepFitAffine(dat, model, opt)

    % Detect parallelisation scheme
    % ------------------------
    if strcmpi(opt.loop, 'subject')
        opt.loop = '';
        opt.par  = 0;
        opt.verbose = false;
    end
        
    if isfield(dat, 'model')
        noisemodel = dat.model;
    else
        noisemodel = opt.model;
    end
    
    % Compute phi/jac (needed for Affine fitting)
    % ---------------
    if isfield(dat, 'v')
        [dat.iphi, dat.phi, dat.jac] = exponentiateVelocity(dat.v, ...
            'iphi', 'phi', 'jac', ...
            'itgr', opt.itgr, 'vs', opt.vs, ...
            'prm', opt.prm, 'debug', opt.debug, ...
            'output', {dat.iphi, dat.phi, dat.jac});
        iphi = dat.iphi;
        phi  = dat.phi;
        jac  = dat.jac;
    else
        iphi = warps('identity', opt.lat);
        phi  = [];
        jac  = [];
    end
    
    % Gauss-Newton iterations
    % -----------------------
    % It is useful to actually find a mode of the posterior (and not only
    % an improved value) when we use the Laplace precision for the update
    % of W. In that case, setting gnit > 1 might help converge faster.
    dat.okq = false;
    for i=1:opt.gnit

        % Compute gradient/hessian
        % ------------------------
        
        [dat.gq, dat.hq] = ghMatchingAffine(noisemodel, ...
            model.mu, dat.pf, dat.c, ...
            model.gmu, dat.A, opt.affine_basis, ...
            phi, jac, 'bb', dat.bb, ...
            'Mmu', model.Mmu, 'loop', opt.loop, 'par', opt.par, ...
            'debug', opt.debug, 'approx', opt.happrox);
        
        if checkarray(model.regq)
            rind = opt.affine_rind;
            [gq, hq] = ghPriorAffine(dat.q(rind), model.regq, 'debug', opt.debug);
            dat.gq(rind)      = dat.gq(rind)      + gq;
            dat.hq(rind,rind) = dat.hq(rind,rind) + hq;
            clear gq hq
        else
            rind = [];
        end
        
        dat.hq = loadDiag(dat.hq); % Additional regularisation for robustness)

        % Compute search direction
        % ------------------------
        dq = -dat.hq \ dat.gq;

        % Line search
        % -----------
%         dat = oneUpdate(dat, model, opt, struct('llmp', true));
        [okq, q, llm, ~, A, pf, c, bb] = lsAffine(...
            noisemodel, dq, dat.q, dat.llm, model.mu, dat.f, ...
            'B', opt.affine_basis, 'regq', model.regq, 'rind', rind, ...
            'iphi', iphi, ...
            'Mf', dat.Mf, 'Mmu', model.Mmu, 'nit', opt.lsit, ...
            'par', opt.par, 'verbose', opt.verbose, 'debug', opt.debug);

        % Store better values
        % -------------------
        dat.okq = dat.okq || okq;
        if okq
            dat.q       = q;
            dat.qq      = q * q';
            dat.llm     = llm;
            dat.A       = A;
            dat.pf      = prepareOnDisk(dat.pf, size(pf));
            dat.pf(:)   = pf(:);
            dat.c       = prepareOnDisk(dat.c, size(c));
            dat.c(:)    = c(:);
            dat.bb      = bb;
%             dat.ipsi    = prepareOnDisk(dat.ipsi, size(ipsi));
%             dat.ipsi(:) = ipsi(:);
        else
            break
        end

%         dat = oneUpdate(dat, model, opt, ...
%                         struct('wmu', true, 'llmw', true), ...
%                         struct('wmu', true));
        
    end
        
    % Compute Laplace covariance
    % --------------------------

    if okq
        dat.hq = ghMatchingAffine(noisemodel, ...
            model.mu, dat.pf, dat.c, ...
            model.gmu, dat.A, opt.affine_basis, ...
            dat.phi, dat.jac, ...
            'bb', dat.bb', 'hessian', true, ...
            'Mmu', model.Mmu, 'loop', opt.loop, 'par', opt.par, ...
            'debug', opt.debug, 'approx', opt.happrox);

        if checkarray(model.regq)
            [~, hq] = ghPriorAffine(dat.q(rind), model.regq, 'debug', opt.debug);
            dat.hq(rind,rind) = dat.hq(rind,rind) + hq;
            clear hq
        end

        dat.hq = loadDiag(dat.hq); % Additional regularisation for robustness
        
    end
    dat.Sq = inv(dat.hq);
    
    % Cleaning
    % --------
    % I should probably clear variables and remove files that are not
    % useful anymore. This will cause less disk and broadband usage.
    toclean = {'gq', 'hq', 'iphi', 'ipsi', 'phi', 'jac'};
    for i=1:numel(toclean)
        field = toclean{i};
        dat.(field) = rmarray(dat.(field));
    end
end

function [dat, model] = batchFitAffine(dat, model, opt)

    % --- Detect parallelisation scheme
    if strcmpi(opt.loop, 'subject') && opt.par > 0
        batch = opt.batch;
    else
        batch = 1;
    end

    % Init gradient/hessian
    % ---------------------
    model.llm  = 0;
    okq        = 0;
    M          = size(opt.affine_basis, 3);
    model.q    = zeros(M, 1);
    model.qq   = zeros(M, M);
    model.Sq   = zeros(M, M);
    
    % --- Batch processing
    if opt.verbose, before = plotBatchBegin('Fit Affine'); end;
    for i=1:ceil(opt.N/batch)
        n1 = (i-1)*batch + 1;
        ne = min(opt.N, i*batch);

        if opt.verbose, before = plotBatch(i, batch, opt.N, 50, before); end;
    
        % Compute subjects grad/hess w.r.t. initial velocity
        % --------------------------------------------------
        dat(n1:ne) = distribute('oneStepFitAffine', dat(n1:ne), model, opt);
        
        for n=n1:ne
            
            % Add individual contributions
            % ----------------------------
            model.llm  = model.llm + dat(n).llm;
            okq        = okq       + dat(n).okq;
            model.q    = model.q   + dat(n).q;
            model.qq   = model.qq  + dat(n).qq;
            model.Sq   = model.Sq  + dat(n).Sq;
        end
        
    end
    if opt.verbose, fprintf(' | %4d / %4d', okq, opt.N); end;
    if opt.verbose, plotBatchEnd; end;

end

% ------------------
% Latent coordinates
% ------------------

function dat = oneStepFitLatent(dat, model, opt)

    % Detect parallelisation scheme
    % ------------------------
    if strcmpi(opt.loop, 'subject')
        opt.loop = '';
        opt.par  = 0;
        opt.verbose = false;
    end
    
    if isfield(dat, 'model')
        noisemodel = dat.model;
    else
        noisemodel = opt.model;
    end
    
    % Gauss-Newton iterations
    % -----------------------
    % It is useful to actually find a mode of the posterior (and not only
    % an improved value) when we use the Laplace precision for the update
    % of W. In that case, setting gnit > 1 might help converge faster.
    dat.okz = false;
    for i=1:opt.gnit

        % Compute gradient/hessian
        % ------------------------

        [dat.gz, dat.hz] = ghMatchingLatent(noisemodel, ...
            model.mu, dat.pf, dat.c, model.gmu, model.w, ...
            'bb', dat.bb, ...
            'loop', opt.loop, 'par', opt.par, 'debug', opt.debug, ...
            'output', {dat.gz, dat.hz});

        [gz, hz] = ghPriorLatent(dat.z, model.regz, 'debug', opt.debug);
        dat.gz = dat.gz + gz;
        dat.hz = dat.hz + hz;
        clear gz hz

        dat.hz = loadDiag(dat.hz); % Additional regularisation for robustness

        % Compute search direction
        % ------------------------
        dz = -dat.hz \ dat.gz;

        % Line search
        % -----------
        if isfield(dat, 'A'),         A = dat.A;
        else                          A = eye(4); end
%         dat = oneUpdate(dat, model, opt, struct('llmp', true));
        [okz, z, llm, ~, v, ~, pf, c, bb] = lsLatent(...
            noisemodel, dz, dat.z, dat.v, dat.llm, ...
            model.w, model.mu, dat.f, ...
            'regz', model.regz, ...
            'A', A, 'Mf', dat.Mf, 'Mmu', model.Mmu, ...
            'nit', opt.lsit, 'itgr', opt.itgr, 'prm', opt.prm);

        % Store better values
        % -------------------
        dat.okz = dat.okz || okz;
        if okz
            dat.z       = z;
            dat.zz      = z * z';
            dat.llm     = llm;
            dat.v(:)    = v(:);
%             dat.iphi(:) = iphi(:);
            dat.pf      = prepareOnDisk(dat.pf, size(pf));
            dat.pf(:)   = pf(:);
            dat.c       = prepareOnDisk(dat.c, size(c));
            dat.c(:)    = c(:);
            dat.bb      = bb;
%             dat.ipsi(:) = ipsi(:);
        else
            break
        end
        
%         dat = oneUpdate(dat, model, opt, ...
%                         struct('wmu', true, 'llmw', true), ...
%                         struct('wmu', true));
    end
    
    % Compute Laplace covariance
    % --------------------------

    if okz
        dat.hz = ghMatchingLatent(noisemodel, ...
            model.mu, dat.pf, dat.c, model.gmu, model.w,...
            'bb', dat.bb, 'hessian', true);

        [~, hz] = ghPriorLatent(dat.z, model.regz);
        dat.hz = dat.hz + hz;
        clear hz

        dat.hz = loadDiag(dat.hz); % Additional regularisation for robustness
        
    end
    dat.Sz = inv(dat.hz);
        
    % Cleaning
    % --------
    % I should probably clear variables and remove files that are not
    % useful anymore. This will cause less disk and broadband usage.
    toclean = {'gz', 'hz', 'iphi', 'ipsi'};
    for i=1:numel(toclean)
        field = toclean{i};
        dat.(field) = rmarray(dat.(field));
    end
end

function [dat, model] = batchFitLatent(dat, model, opt)

    % --- Detect parallelisation scheme
    if strcmpi(opt.loop, 'subject') && opt.par > 0
        batch = opt.batch;
    else
        batch = 1;
    end

    % Init gradient/hessian
    % ---------------------
    model.llm  = 0;
    model.z    = zeros(opt.K, 1);
    model.zz   = zeros(opt.K, opt.K);
    model.Sz   = zeros(opt.K, opt.K);
    okz        = 0;
    
    % --- Batch processing
    if opt.verbose, before = plotBatchBegin('Fit Latent'); end;
    for i=1:ceil(opt.N/batch)
        n1 = (i-1)*batch + 1;
        ne = min(opt.N, i*batch);

        if opt.verbose, before = plotBatch(i, batch, opt.N, 50, before); end;
    
        % Compute subjects grad/hess w.r.t. initial velocity
        % --------------------------------------------------
        dat(n1:ne) = distribute('oneStepFitLatent', dat(n1:ne), model, opt);
        
        for n=n1:ne
            
            % Add individual contributions
            % ----------------------------
            model.z    = model.z   + dat(n).z;
            model.zz   = model.zz  + dat(n).zz;
            model.Sz   = model.Sz  + dat(n).Sz;
            model.llm  = model.llm + dat(n).llm;
            okz        = okz + dat(n).okz;
            
        end
        
    end
    if opt.verbose, fprintf(' | %4d / %4d', okz, opt.N); end;
    if opt.verbose, plotBatchEnd; end;
end

% --------------
% Residual field
% --------------

function dat = oneStepFitResidual(dat, model, opt)

    % Detect parallelisation scheme
    % ------------------------
    if strcmpi(opt.loop, 'subject')
        opt.loop = '';
        opt.par  = 0;
        opt.verbose = false;
    end
    
    if isfield(dat, 'model')
        noisemodel = dat.model;
    else
        noisemodel = opt.model;
    end
    
    if isfield(dat, 'A'),         A = dat.A;
    else                          A = eye(4); end
    if isfield(model, 'lambda'),  lambda = model.lambda;
    else                          lambda = 1; end
    
    % Gauss-Newton iterations
    % -----------------------
    % It is useful to actually find a mode of the posterior (and not only
    % an improved value) when we use the Laplace precision for the update
    % of W. In that case, setting gnit > 1 might help converge faster.
    dat.okr = false;
    for i=1:opt.gnit

        % Compute gradient/hessian
        % ------------------------
        lat = [size(model.mu) 1];
        dat.gr = prepareOnDisk(dat.gr, [lat(1:3) 3]); 
        dat.gr(:) = 0;
        dat.hr = prepareOnDisk(dat.hr, [lat(1:3) 6]); 
        dat.hr(:) = 0;
        bx = dat.bb.x;
        by = dat.bb.y;
        bz = dat.bb.z;
        
        [dat.gr(bx,by,bz,:), dat.hr(bx,by,bz,:)] = ghMatchingVel(...
            noisemodel, ...
            model.mu, dat.pf, dat.c, model.gmu, 'bb', dat.bb, ...
            'loop', opt.loop, 'par', opt.par, 'debug', opt.debug);
        
        dat.gr(:,:,:,:) = dat.gr(:,:,:,:) + ghPriorVel(dat.r, ...
            sqrt(sum(model.Mmu(1:3,1:3).^2)), lambda * opt.prm);

        % Compute search direction
        % ------------------------
        dr = -spm_diffeo('fmg', single(dat.hr), single(dat.gr), ...
            double([sqrt(sum(model.Mmu(1:3,1:3).^2)) lambda * opt.prm 2 2]));

        % Line search
        % -----------
%         dat = oneUpdate(dat, model, opt, struct( 'llmp', true));
        [okr, r, llm, llr, ~, pf, c, bb, ~, v] = lsVelocity(...
            noisemodel, dr, dat.r, dat.llm, ...
            model.mu, dat.f, 'v0', dat.v, 'A', A, ....
            'Mf', dat.Mf, 'Mmu', model.Mmu, 'nit', opt.lsit, ...
            'itgr', opt.itgr, 'prm', lambda * opt.prm, ...
            'par', opt.par, 'verbose', opt.verbose, 'debug', opt.debug);

        % Store better values
        % -------------------
        dat.okr = dat.okr || okr;
        if okr
            dat.r(:)    = r(:);
            dat.llm     = llm;
            dat.llr     = llr;
            dat.v(:)    = v(:);
%             dat.iphi(:) = iphi(:);
            dat.pf      = prepareOnDisk(dat.pf, size(pf));
            dat.pf(:)   = pf(:);
            dat.c       = prepareOnDisk(dat.c, size(c));
            dat.c(:)    = c(:);
            dat.bb      = bb;
%             dat.ipsi(:) = ipsi(:);
        else
            break
        end
        
%         dat = oneUpdate(dat, model, opt, ...
%                         struct('wmu', true, 'llmw', true), ...
%                         struct('wmu', true));
    end
    
    % Compute Laplace covariance
    % --------------------------

    if okr
        bx = dat.bb.x;
        by = dat.bb.y;
        bz = dat.bb.z;
        dat.hr(bx,by,bz,:) = ghMatchingVel(noisemodel, ...
            model.mu, dat.pf, dat.c, model.gmu, ...
            'bb', dat.bb, 'hessian', true, ...
            'loop', opt.loop, 'par', opt.par, 'debug', opt.debug);
    end
    
    % Statistics for precision update
    % -------------------------------
    
    % err: trace(E[RR]L) -> for lambda upadte
    % klr1: prior part of the KL divergence
    % klr2: posterior part of the KL divergence
    
    lam = model.lambda;
    hr = single(numeric(dat.hr));
    r  = single(numeric(dat.r));
    vs = sqrt(sum(model.Mmu(1:3,1:3).^2));
    
    % 1) Compute all elements
    
    % - trace
    tr = spm_diffeo('trapprox', single(hr/lam), double([vs opt.prm]));
    tr = tr(1);
    % - reg prior
    llr = llPriorVelocity(r, 'fast', 'vs', vs,  'prm', opt.prm);
    % - det prior
    [~, ld1] = spm_shoot_greens('kernel', double(opt.lat), double([vs opt.prm]));
    ld1 = ld1(1);
    % - det posterior
    K = spm_diffeo('kernel', double(opt.lat), double([vs opt.prm]));
    hr(:,:,:,1) = hr(:,:,:,1) + lam*K(1,1,1,1,1);
    hr(:,:,:,2) = hr(:,:,:,2) + lam*K(1,1,1,2,2);
    if size(hr, 3) == 1
        hr(:,:,:,3)   = 1;
        hr(:,:,:,5:6) = 0;
    else
        hr(:,:,:,3) = hr(:,:,:,3) + lam*K(1,1,1,3,3);
    end
    ld2 = sumall(log(abs(pointwise3(hr, 'd'))));
    
    % 2) Sum each statistic
    dat.err  = 0.5*(tr/lam + llr);
    dat.trr  = tr;   % (keep track so that we do not need recomputing it)
    dat.llr  = llr;  % (keep track so that we do not need recomputing it)
    dat.klr1 = 0.5*(tr - prod(opt.lat)*3*log(lam) + lam*llr);
    dat.klr2 = 0.5*(ld2 - ld1 - prod(opt.lat)*3);
    dat.klr  = dat.klr1 + dat.klr2;
    
    % Cleaning
    % --------
    % I should probably clear variables and remove files that are not
    % useful anymore. This will cause less disk and broadband usage.
    toclean = {'gr', 'hr', 'iphi', 'ipsi'};
    for i=1:numel(toclean)
        field = toclean{i};
        dat.(field) = rmarray(dat.(field));
    end
end

function [dat, model] = batchFitResidual(dat, model, opt)

    % --- Detect parallelisation scheme
    if strcmpi(opt.loop, 'subject') && opt.par > 0
        batch = opt.batch;
    else
        batch = 1;
    end

    % Init gradient/hessian
    % ---------------------
    model.llm  = 0;
    model.llr  = 0;
    model.err  = 0;
    okr        = 0;
    
    % --- Batch processing
    if opt.verbose, before = plotBatchBegin('Fit Res'); end;
    for i=1:ceil(opt.N/batch)
        n1 = (i-1)*batch + 1;
        ne = min(opt.N, i*batch);

        if opt.verbose, before = plotBatch(i, batch, opt.N, 50, before); end;
    
        % Compute subjects grad/hess w.r.t. initial velocity
        % --------------------------------------------------
        dat(n1:ne) = distribute('oneStepFitResidual', dat(n1:ne), model, opt);
        
        
        for n=n1:ne
            
            % Add individual contributions
            % ----------------------------
            model.llm  = model.llm + dat(n).llm;
            model.llr  = model.llr + dat(n).llr;
            model.err  = model.err + dat(n).err;
            okr        = okr + dat(n).okr;
            
        end
        
    end
    if opt.verbose, fprintf(' | %4d / %4d', okr, opt.N); end;
    if opt.verbose, plotBatchEnd; end;

end

% -------------------------------------------------------------------------
%    M-step
% -------------------------------------------------------------------------

function dat = oneStepGradHessVelocity(dat, model, opt)

    if isfield(dat, 'model')
        noisemodel = dat.model;
    else
        noisemodel = opt.model;
    end
    
    [dat.gv, dat.hv] = ghMatchingVel(noisemodel, ...
        model.mu, dat.pf, dat.c, model.gmu, ...
        'bb', dat.bb, 'output', {dat.gv, dat.hv});
end

% --------
% Subspace
% --------

function [dat, model] = batchGradHessSubspace(dat, model, opt)

    % --- Detect parallelisation scheme
    if strcmpi(opt.loop, 'subject') && opt.par > 0
        batch = opt.batch;
    else
        batch = 1;
    end

    % Init gradient/hessian
    % ---------------------
    dim     = [size(model.w) 1 1 1];
    model.gw = prepareOnDisk(model.gw, [dim(1:3) 3 dim(5)], 'type', 'float32');
    model.hw = prepareOnDisk(model.hw, [dim(1:3) 6 dim(5)], 'type', 'float32');
    for k=1:dim(5)
        model.gw(:,:,:,:,k) = 0;
        model.hw(:,:,:,:,k) = 0;
    end
    
    % --- Batch processing
    if opt.verbose, before = plotBatchBegin('GH PG'); end;
    for i=1:ceil(opt.N/batch)
        n1 = (i-1)*batch + 1;
        ne = min(opt.N, i*batch);

        if opt.verbose, before = plotBatch(i, batch, opt.N, 50, before); end;
    
        % Compute subjects grad/hess w.r.t. initial velocity
        % --------------------------------------------------
        dat(n1:ne) = distribute('oneStepGradHessVelocity', dat(n1:ne), model, opt);
        
        for n=n1:ne
            
            % Add individual contributions
            % ----------------------------
            for k=1:opt.K
                bx = dat(n).bb.x;
                by = dat(n).bb.y;
                bz = dat(n).bb.z;
                model.gw(bx,by,bz,:,k) = model.gw(bx,by,bz,:,k) + numeric(dat(n).gv) * single(dat(n).z(k));
                model.hw(bx,by,bz,:,k) = model.hw(bx,by,bz,:,k) + numeric(dat(n).hv) * single(dat(n).z(k))^2;
            end
            
            % Clear individual grad/hess
            % --------------------------
            dat(n).gv = rmarray(dat(n).gv);
            dat(n).hv = rmarray(dat(n).hv);
        end
        
    end
    if opt.verbose, plotBatchEnd; end;

end

% -------------------------------------------------------------------------
%    Update
% -------------------------------------------------------------------------

function dat = oneUpdate(dat, model, opt, todo, toclean)

    if nargin < 5
        toclean = struct;
        if nargin < 4
            todo = struct;
        end
    end

    if strcmpi(opt.loop, 'subject')
        loop = '';
        par  = 0;
    else
        loop = opt.loop;
        par  = opt.par;
    end
    
    if isfield(dat, 'model')
        noisemodel = dat.model;
    else
        noisemodel = opt.model;
    end
    
    % --- Transform
    if isfield(todo, 'v')
        if isfield(dat, 'z') && isfield(dat, 'r')
            dat.v = reconstructVelocity('latent', dat.z, 'subspace', model.w, ...
                'residual', dat.r, ...
                'debug', opt.debug, 'output', dat.v, ...
                'loop', loop, 'par', par);
        elseif isfield(dat, 'z')
            dat.v = reconstructVelocity('latent', dat.z, 'subspace', model.w, ...
                'debug', opt.debug, 'output', dat.v, ...
                'loop', loop, 'par', par);
        elseif isfield(dat, 'r')
            dat.v = reconstructVelocity('residual', dat.r, ...
                'debug', opt.debug, 'output', dat.v, ...
                'loop', loop, 'par', par);
        else
            error('Need at least latent space or residual field to reconstruct velocity')
        end
    end
    if isfield(todo, 'iphi') && isfield(todo, 'phi') && isfield(todo, 'jac')
        [dat.iphi, dat.phi, dat.jac] = exponentiateVelocity(dat.v, ...
            'iphi', 'phi', 'jac', ...
            'itgr', opt.itgr, 'vs', opt.vs, ...
            'prm', opt.prm, 'debug', opt.debug, ...
            'output', {dat.iphi, dat.phi, dat.jac});
    elseif isfield(todo, 'iphi')
        dat.iphi = exponentiateVelocity(dat.v, 'iphi', ...
            'itgr', opt.itgr, 'vs', opt.vs, ...
            'prm', opt.prm, 'debug', opt.debug, 'output', dat.iphi);
    elseif isfield(todo, 'phi') && isfield(todo, 'jac')
        [dat.phi, dat.jac] = exponentiateVelocity(dat.v, ...
            'phi', 'jac', ...
            'itgr', opt.itgr, 'vs', opt.vs, ...
            'prm', opt.prm, 'debug', opt.debug, ...
            'output', {dat.phi, dat.jac});
    end
    if isfield(toclean, 'v')
        dat.v = rmarray(dat.v);
    end
    if isfield(todo, 'A')
        dat.A = exponentiateAffine(dat.q, opt.affine_basis, ...
            'debug', opt.debug, 'output', dat.A);
    end
    if isfield(todo, 'ipsi')
        latf = [size(dat.f) 1];
        latf = latf(1:3);
        if isfield(dat, 'A') && isfield(dat, 'iphi')
            dat.ipsi = reconstructIPsi(dat.A, dat.iphi, ...
                'lat', latf, 'Mf', dat.Mf, 'Mmu', model.Mmu, ...
                'output', dat.ipsi, 'debug', opt.debug);
        elseif isfield(dat, 'iphi')
            dat.ipsi = reconstructIPsi(eye(4), dat.iphi, ...
                'lat', latf, 'Mf', dat.Mf, 'Mmu', model.Mmu, ...
                'output', dat.ipsi, 'debug', opt.debug);
        elseif isfield(dat, 'A')
            dat.ipsi = reconstructIPsi(dat.A, warps('identity', opt.lat), ...
                'lat', latf, 'Mf', dat.Mf, 'Mmu', model.Mmu, ...
                'output', dat.ipsi, 'debug', opt.debug);
            error('Only affine is not implemented yet')
        else
            error('Need at least iphi or A to reconstruct transformation')
        end
    end
    if isfield(toclean, 'iphi')
        dat.iphi = rmarray(dat.iphi);
    end
    
    % --- Push/Warp
    if isfield(todo, 'pf') || isfield(todo, 'c')
        [dat.pf, dat.c, dat.bb] = pushImage(dat.ipsi, dat.f, opt.lat, ...
            'loop', loop, 'par', par, ...
            'output', {dat.pf, dat.c}, 'debug', opt.debug);
    end
    if isfield(todo, 'wmu')
        if opt.tpm
            a = warp(dat.ipsi, model.a, opt.itrp, opt.bnd, ...
                'par', par, 'debug', opt.debug);
            dat.wmu = reconstructProbaTemplate(a, ...
                'loop', loop, 'par', par, 'debug', opt.debug, ...
                'output', dat.wmu);
        else
            dat.wmu = warp(dat.ipsi, model.mu, opt.itrp, opt.bnd, ...
                'output', dat.wmu, 'debug', opt.debug);
        end
    end
    if isfield(toclean, 'ipsi')
        dat.ipsi = rmarray(dat.ipsi);
    end
    if isfield(todo, 'llm') || isfield(todo, 'llmp')
        % Pushed image
        dat.llm = llMatching(noisemodel, model.mu, dat.pf, dat.c, ...
            'bb', dat.bb, 'loop', loop, 'par', par, ...
            'debug', opt.debug);
    elseif isfield(todo, 'llmw')
        % Warped template
        dat.llm = llMatching(noisemodel, dat.wmu, dat.f, ...
            'loop', loop, 'par', par, ...
            'debug', opt.debug);
    end
    if isfield(toclean, 'wmu')
        dat.wmu = rmarray(dat.wmu);
    end
    
    % --- G/H affine
    if isfield(todo, 'gq') && isfield(todo, 'hq')
        [dat.gq, dat.hq] = ghMatchingAffine(noisemodel, ...
            model.mu, dat.pf, dat.c, model.gmu, dat.A, ...
            opt.affine_basis, dat.phi, dat.jac, ...
            'bb', dat.bb, 'loop', loop, 'par', par);

        [gq, hq] = ghPriorAffine(dat.q, model.regq);
        dat.gq = dat.gq + gq;
        dat.hq = dat.hq + hq;
        clear gq hq
    elseif isfield(todo, 'hq')
        dat.hq = ghMatchingAffine(noisemodel, ...
            model.mu, dat.pf, dat.c, model.gmu, dat.A, ...
            opt.affine_basis, dat.phi, dat.jac, ...
            'bb', dat.bb', 'hessian', true, 'loop', loop, 'par', par);

        [~, hq] = ghPriorAffine(dat.q, model.regq);
        dat.hq = dat.hq + hq;
        clear hq
    elseif isfield(todo, 'gq')
        dat.gq = ghMatchingAffine(noisemodel, ...
            model.mu, dat.pf, dat.c, model.gmu, dat.A, ...
            opt.affine_basis, dat.phi, dat.jac, ...
            'bb', dat.bb, 'loop', loop, 'par', par);

        gq = ghPriorAffine(dat.q, model.regq);
        dat.gq = dat.gq + gq;
        clear gq
    end
    if isfield(todo, 'hq')
        dat.hq = loadDiag(dat.hq);
    end
    if isfield(toclean, 'A')
        dat.A = rmarray(dat.A);
    end
    
    % --- LL affine
    if isfield(todo, 'lllq')
        dat.lllq = llLaplace(dat.hq);
    end
    if isfield(todo, 'Sq')
        dat.Sq = inv(dat.hq);
    end
    if isfield(toclean, 'gq')
        dat.gq = rmarray(dat.gq);
    end
    if isfield(toclean, 'hq')
        dat.hq = rmarray(dat.hq);
    end
    if isfield(todo, 'llq')
        dat.llq = llPriorLatent(dat.q, model.regq, 'debug', opt.debug);
    end
    
    % --- G/H latent
    if isfield(todo, 'gz') && isfield(todo, 'hz')
        [dat.gz, dat.hz] = ghMatchingLatent(noisemodel, ...
            model.mu, dat.pf, dat.c, model.gmu, model.w, ...
            'bb', dat.bb, 'loop', loop, 'par', par);

        [gz, hz] = ghPriorLatent(dat.z, model.regz);
        dat.gz = dat.gz + gz;
        dat.hz = dat.hz + hz;
        clear gz hz
    elseif isfield(todo, 'hz')
        dat.hz = ghMatchingLatent(noisemodel, ...
            model.mu, dat.pf, dat.c, model.gmu, model.w, ...
            'bb', dat.bb, 'hessian', true, 'loop', loop, 'par', par);

        [~, hz] = ghPriorLatent(dat.z, model.regz);
        dat.hz = dat.hz + hz;
        clear hz
    elseif isfield(todo, 'gz')
        dat.gz = ghMatchingLatent(noisemodel, ...
            model.mu, dat.pf, dat.c, model.gmu, model.w, ...
            'bb', dat.bb, 'loop', loop, 'par', par);

        gz = ghPriorLatent(dat.z, model.regz);
        dat.gz = dat.gz + gz;
        clear gz
    end
    if isfield(todo, 'hz')
        dat.hz = loadDiag(dat.hz);
    end
    
    % --- LL latent
    if isfield(todo, 'lllz')
        dat.lllz = llLaplace(dat.hz);
    end
    if isfield(todo, 'Sz')
        dat.Sz = inv(dat.hz);
    end
    if isfield(toclean, 'gz')
        dat.gz = rmarray(dat.gz);
    end
    if isfield(toclean, 'hz')
        dat.hz = rmarray(dat.hz);
    end
    if isfield(todo, 'llz')
        dat.llz = llPriorLatent(dat.z, model.regz, 'debug', opt.debug);
    end
    
    % --- G/H residual
    if isfield(todo, 'hr')
        dat.hr = prepareOnDisk(dat.hr, [opt.lat(1:3) 6]);
        dat.hr(:) = 0;
        bx = dat.bb.x;
        by = dat.bb.y;
        bz = dat.bb.z;
        dat.hr(bx,by,bz,:) = ghMatchingVel(noisemodel, ...
            model.mu, dat.pf, dat.c, model.gmu, ...
            'bb', dat.bb, 'hessian', true, ...
            'loop', loop, 'par', par, 'debug', opt.debug);
    end
    if isfield(todo, 'klr')
        lam  = model.lambda;
        lam0 = model.lambda_prev;
        hr   = single(numeric(dat.hr));
        r    = single(numeric(dat.r));
        vs   = sqrt(sum(model.Mmu(1:3,1:3).^2));

        % 1) Compute all elements

        % - trace
        tr = spm_diffeo('trapprox', single(hr/lam0), double([vs opt.prm]));
        tr = tr(1);
        % - reg prior
        llr = llPriorVelocity(r, 'fast', 'vs', vs,  'prm', opt.prm);
        clear r
        % - det prior
        [~, ld1] = spm_shoot_greens('kernel', double(opt.lat), double([vs opt.prm]));
        ld1 = ld1(1);
        % - det posterior
        K = spm_diffeo('kernel', double(opt.lat), double([vs opt.prm]));
        hr(:,:,:,1) = hr(:,:,:,1) + lam*K(1,1,1,1,1);
        hr(:,:,:,2) = hr(:,:,:,2) + lam*K(1,1,1,2,2);
        if size(hr, 3) == 1
            hr(:,:,:,3)   = 1;
            hr(:,:,:,5:6) = 0;
        else
            hr(:,:,:,3) = hr(:,:,:,3) + lam*K(1,1,1,3,3);
        end
        ld2 = sumall(log(abs(pointwise3(hr, 'd'))));
        clear hr

        % 2) Sum each statistic
        dat.trr  = tr;   % (keep track so that we do not need recomputing it)
        dat.llr  = llr;  % (keep track so that we do not need recomputing it)
        dat.klr1 = 0.5*(lam/lam0*tr - prod(opt.lat)*3*log(lam) + lam*llr);
        dat.klr2 = 0.5*(ld2 - ld1 - prod(opt.lat)*3);
        dat.klr  = dat.klr1 + dat.klr2;
    end
    if isfield(todo, 'llr')
        dat.llr = llPriorVelocity(dat.r, ...
            'vs', sqrt(sum(model.Mmu(1:3,1:3).^2)), ...
            'prm', opt.prm);
    end
    if isfield(todo, 'klru')
        lam  = model.lambda;
        lam0 = model.lambda_prev;
        dat.klr1 = 0.5*(lam/lam0*dat.trr - prod(opt.lat)*3*log(lam) + lam*dat.llr);
        dat.klr  = dat.klr1 + dat.klr2;
    end
    if isfield(toclean, 'hr')
        dat.hr = rmarray(dat.hr);
    end
    
    % --- Clean pushed
    if isfield(toclean, 'pf')
        dat.pf = rmarray(dat.pf);
    end
    if isfield(toclean, 'c')
        dat.c = rmarray(dat.c);
    end
end

function dat = batchUpdate(dat, model, opt, todo, varargin)

    % --- Parse inputs
    p = inputParser;
    p.FunctionName = 'batchUpdate';
    p.addRequired('dat',    @isstruct);
    p.addRequired('model',  @isstruct);
    p.addRequired('opt',    @isstruct);
    p.addRequired('todo',   @(X) ischar(X) || iscell(X));
    p.addParameter('clean',  {}, @(X) ischar(X) || iscell(X));
    p.parse(dat, model, opt, todo, varargin{:});
    clean = p.Results.clean;
    
    if opt.debug, fprintf('* batchUpdate\n'); end;
    
    % --- Detect parallelisation scheme
    if strcmpi(opt.loop, 'subject') && opt.par > 0
        batch = opt.batch;
    else
        batch = 1;
    end
    
    % --- Transform todo into a structure for simplicity
    if ~iscell(todo)
        todo = {todo};
    end
    stodo = struct;
    for i=1:numel(todo)
        stodo.(todo{i}) = true;
    end
    if ~iscell(clean)
        clean = {clean};
    end
    stoclean = struct;
    for i=1:numel(clean)
        stoclean.(clean{i}) = true;
    end
    
    % --- Batch processing
    if opt.verbose, before = plotBatchBegin('Update'); end;
    for i=1:ceil(numel(dat)/batch)
        n1 = (i-1)*batch + 1;
        ne = min(numel(dat), i*batch);
        
        if opt.verbose, before = plotBatch(i, batch, numel(dat), 50, before); end;
        
        dat(n1:ne) = distribute('OneUpdate', dat(n1:ne), model, opt, stodo, stoclean);
    end
    if opt.verbose, plotBatchEnd; end;
    
end

% =========================================================================

function dat = distribute(func, dat, model, opt, varargin)

    % --- Convert function name to function handle
    if ischar(func)
        funcname = func;
        func = @(varargin) batchProcess(funcname, varargin{:});
        funcstr = ['@(varargin) batchProcess(''' funcname ''', varargin{:})']; 
    else
        funcstr = func2str(func);
    end

    % --- No distribution
    if ~strcmpi(opt.loop, 'subject')
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
        opt.loop    = '';
        opt.par     = false;
        opt.verbose = true;
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
            if opt.par
                parfor (n=1:numel(dat), double(opt.par))
                    dat(n) = func(dat(n), model, opt, varargin{:});
                end
            else
                for n=1:numel(dat)
                    dat(n) = func(dat(n), model, opt, varargin{:});
                end
            end
        else
            if opt.par
                parfor (n=1:numel(dat), double(opt.par))
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
        