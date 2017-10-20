function varargout = pgraBatchProcess(id, varargin)
% Collection of functions that require updating some subject-specific
% variables. Depending on the needs, as much as possible is done in order
% to minimise i/o, memory and computational costs.
%
% FORMAT [dat, model] = batchProcess('InitAffine', mode, dat, model, opt)
%
% FORMAT [dat, model] = batchProcess('InitResidual', mode, dat, model, opt)
%
% FORMAT [dat, model] = batchProcess('InitLatent', mode, dat, model, opt)
% > Initialise latent coordinates and insure they are zero-centered (but
%   orthogonalisation is not done here).
%   Update>> dat.z, dat.zz, dat.Sz, model.z, model.zz, model.Sz
%
% FORMAT dat = batchProcess('RotateLatent', dat, opt, R)
% > Apply a rotation matrix to all latent coordinates
%   Update>> dat.z, dat.zz, dat.S
%
% FORMAT dat = batchProcess('CentreLatent', dat, opt, (mean))
% > Subtract the mean to all latent coordinates
%   Update>> dat.z, dat.zz
%
% FORMAT [dat, model] = batchProcess('FitAffine', dat, model, opt)
%
% FORMAT [dat, model] = batchProcess('FitLatent', dat, model, opt)
%
% FORMAT [dat, model] = batchProcess('FitResidual', dat, model, opt)
%
% FORMAT [dat, model] = batchProcess('GradHessSubspace', dat, model, opt)
%
% FORMAT [dat, model] = batchProcess('GradHessSigma', dat, model, opt)
%
% FORMAT dat = batchProcess('Update', dat, model, opt, {...}, 'clean', {...})
% > Generic tool to update specified subject variables. Names of variables
%   to update are provided in the first list, and names of variables to
%   discard in the end (to save memory/disk) can be provided in the second
%   list.
%   The model is never updated, consequently, the processing is not really
%   performed by batch but distributed as specified by opt.loop/opt.par.


    switch lower(id)
        case 'initaffine'
            [varargout{1:nargout}] = batchInitAffine(varargin{:});
        case 'initresidual'
            [varargout{1:nargout}] = batchInitResidual(varargin{:});
        case 'initlatent'
            [varargout{1:nargout}] = batchInitLatent(varargin{:});
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
        case 'gradhesssigma'
            [varargout{1:nargout}] = batchGradHessSigma(varargin{:});
        case 'update'
            [varargout{1:nargout}] = batchUpdate(varargin{:});
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

function [dat, model] = batchInitResidual(mode, dat, model, opt)

    % --- Don't parallelise here, there is no need
    
    switch lower(mode)
        case 'zero'
            create = @zeros;
        case 'rand'
            create = @randn;
        otherwise
            error('Unknwon mode %s', mode)
    end
    
    model.llr = 0;
    
    % Init subjects
    % -------------
    if opt.verbose, before = plotBatchBegin('Init R'); end;
    for n=1:opt.N
        if opt.verbose, before = plotBatch(n, 1, opt.N, 50, before); end;
        r        = create([opt.lat 3]);
        dat(n).r = saveOnDisk(dat(n).r, r);
        dat(n).llr = llPriorVelocity(r, ...
            'vs', sqrt(sum(model.Mmu(1:3,1:3).^2)), ...
            'prm', opt.prm, 'debug', opt.debug);
        model.llr = model.llr + dat(n).llr;
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
    
    % Compute phi/jac (needed for Affine fitting)
    % ---------------
    [dat.iphi, dat.phi, dat.jac] = exponentiateVelocity(dat.v, ...
        'iphi', 'phi', 'jac', ...
        'itgr', opt.itgr, 'vs', opt.vs, ...
        'prm', opt.prm, 'debug', opt.debug, ...
        'output', {dat.iphi, dat.phi, dat.jac});
    
    % Gauss-Newton iterations
    % -----------------------
    % It is useful to actually find a mode of the posterior (and not only
    % an improved value) when we use the Laplace precision for the update
    % of W. In that case, setting gnit > 1 might help converge faster.
    dat.okq = false;
    for i=1:opt.gnit

        % Compute gradient/hessian
        % ------------------------
        
        [dat.gq, dat.hq] = ghMatchingAffine(opt.model, ...
            model.mu, dat.pf, dat.c, ...
            model.gmu, dat.A, opt.affine_basis, ...
            dat.phi, dat.jac, ...
            'Mmu', model.Mmu, 'loop', opt.loop, 'par', opt.par, ...
            'debug', opt.debug, 'approx', opt.happrox);
        
        if checkarray(model.regq)
            [gq, hq] = ghPriorAffine(dat.q, model.regq, 'debug', opt.debug);
            dat.gq = dat.gq + gq;
            dat.hq = dat.hq + hq;
            clear gq hq
        end
        
        dat.hq = loadDiag(dat.hq); % Additional regularisation for robustness)

        % Compute search direction
        % ------------------------
        dq = -dat.hq \ dat.gq;

        % Line search
        % -----------
        [okq, q, llm, ~, A, pf, c, ipsi] = lsAffine(...
            opt.model, dq, dat.q, dat.llm, model.mu, dat.f, ...
            'B', opt.affine_basis, 'regq', model.regq, 'iphi', dat.iphi, ...
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
            dat.pf(:)   = pf(:);
            dat.c(:)    = c(:);
            dat.ipsi(:) = ipsi(:);
        else
            break
        end
        
    end
        
    % Compute Laplace covariance
    % --------------------------

    if okq
        [~, dat.hq] = ghMatchingAffine(opt.model, ...
            model.mu, dat.pf, dat.c, ...
            model.gmu, dat.A, opt.affine_basis, ...
            dat.phi, dat.jac, ...
            'Mmu', model.Mmu, 'loop', opt.loop, 'par', opt.par, ...
            'debug', opt.debug, 'approx', opt.happrox);

        if checkarray(model.regq)
            [gq, hq] = ghPriorAffine(dat.q, model.regq, 'debug', opt.debug);
            dat.gq = dat.gq + gq;
            dat.hq = dat.hq + hq;
            clear gq hq
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
        par = opt.par;
        batch = opt.batch;
    else
        par   = 0;
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
        parfor (n=n1:ne, double(par))
            dat(n) = oneStepFitAffine(dat(n), model, opt);
        end
        
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
    
    % Gauss-Newton iterations
    % -----------------------
    % It is useful to actually find a mode of the posterior (and not only
    % an improved value) when we use the Laplace precision for the update
    % of W. In that case, setting gnit > 1 might help converge faster.
    dat.okz = false;
    for i=1:opt.gnit

        % Compute gradient/hessian
        % ------------------------

        [dat.gz, dat.hz] = ghMatchingLatent(opt.model, ...
            model.mu, dat.pf, dat.c, model.gmu, model.w, ...
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
        [okz, z, llm, ~, v, iphi, pf, c, ipsi] = lsLatent(...
            opt.model, dz, dat.z, dat.v, dat.llm, ...
            model.w, model.mu, dat.f, ...
            'regz', model.regz, ...
            'A', dat.A, 'Mf', dat.Mf, 'Mmu', model.Mmu, ...
            'nit', opt.lsit, 'itgr', opt.itgr, 'prm', opt.prm);

        % Store better values
        % -------------------
        dat.okz = dat.okz || okz;
        if okz
            dat.z       = z;
            dat.zz      = z * z';
            dat.llm     = llm;
            dat.v(:)    = v(:);
            dat.iphi(:) = iphi(:);
            dat.pf(:)   = pf(:);
            dat.c(:)    = c(:);
            dat.ipsi(:) = ipsi(:);
        else
            break
        end
        
    end
    
    % Compute Laplace covariance
    % --------------------------

    if okz
        [~, dat.hz] = ghMatchingLatent(opt.model, ...
            model.mu, dat.pf, dat.c, model.gmu, model.w);

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
        par = opt.par;
        batch = opt.batch;
    else
        par   = 0;
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
        parfor (n=n1:ne, double(par))
            dat(n) = oneStepFitLatent(dat(n), model, opt);
        end
        
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
    
    % Gauss-Newton iterations
    % -----------------------
    % It is useful to actually find a mode of the posterior (and not only
    % an improved value) when we use the Laplace precision for the update
    % of W. In that case, setting gnit > 1 might help converge faster.
    dat.okr = false;
    for i=1:opt.gnit

        % Compute gradient/hessian
        % ------------------------
        
        [dat.gr, dat.hr] = ghMatchingVel(opt.model, ...
            model.mu, dat.pf, dat.c, model.gmu, ...
            'loop', opt.loop, 'par', opt.par, 'debug', opt.debug);
        
        dat.gr = dat.gr + ghPriorVel(dat.r, ...
            sqrt(sum(model.Mmu(1:3,1:3).^2)), opt.prm);

        % Compute search direction
        % ------------------------
        dr = -spm_diffeo('fmg', single(dat.hr), single(dat.gr), ...
            [sqrt(sum(model.Mmu(1:3,1:3).^2)) opt.prm 2 2]);

        % Line search
        % -----------
        [okr, r, llm, llr, iphi, pf, c, ipsi, v] = lsVelocity(...
            opt.model, dr, dat.r, dat.llm, ...
            model.mu, dat.f, 'v0', dat.v, 'A', dat.A, ....
            'Mf', dat.Mf, 'Mmu', model.Mmu, 'nit', opt.lsit, ...
            'itgr', opt.itgr, 'prm', opt.prm, 'sigma', model.sigma, ...
            'par', opt.par, 'verbose', opt.verbose, 'debug', opt.debug);

        % Store better values
        % -------------------
        dat.okr = dat.okr || okr;
        if okr
            dat.r       = r;
            dat.llm     = llm;
            dat.llr     = llr;
            dat.v(:)    = v(:);
            dat.iphi(:) = iphi(:);
            dat.pf(:)   = pf(:);
            dat.c(:)    = c(:);
            dat.ipsi(:) = ipsi(:);
        else
            break
        end
        
    end
        
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
        par = opt.par;
        batch = opt.batch;
    else
        par   = 0;
        batch = 1;
    end

    % Init gradient/hessian
    % ---------------------
    model.llm  = 0;
    model.llr  = 0;
    okr        = 0;
    
    % --- Batch processing
    if opt.verbose, before = plotBatchBegin('Fit Res'); end;
    for i=1:ceil(opt.N/batch)
        n1 = (i-1)*batch + 1;
        ne = min(opt.N, i*batch);

        if opt.verbose, before = plotBatch(i, batch, opt.N, 50, before); end;
    
        % Compute subjects grad/hess w.r.t. initial velocity
        % --------------------------------------------------
        parfor (n=n1:ne, double(par))
            dat(n) = oneStepFitResidual(dat(n), model, opt);
        end
        
        for n=n1:ne
            
            % Add individual contributions
            % ----------------------------
            model.llm  = model.llm + dat(n).llm;
            model.llr  = model.llr + dat(n).llr;
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

    [dat.gv, dat.hv] = ghMatchingVel(opt.model, ...
        model.mu, dat.pf, dat.c, model.gmu, ...
        'output', {dat.gv, dat.hv});
end

% --------
% Subspace
% --------

function [dat, model] = batchGradHessSubspace(dat, model, opt)

    % --- Detect parallelisation scheme
    if strcmpi(opt.loop, 'subject') && opt.par > 0
        par = opt.par;
        batch = opt.batch;
    else
        par   = 0;
        batch = 1;
    end

    % Init gradient/hessian
    % ---------------------
    dim     = [size(model.w) 1 1 1];
    model.gw = prepareOnDisk(model.gw, [dim(1:3) 3 dim(5)], 'type', 'float32');
    model.hw = prepareOnDisk(model.hw, [dim(1:3) 6 dim(5)], 'type', 'float32');
    
    % --- Batch processing
    if opt.verbose, before = plotBatchBegin('GH PG'); end;
    for i=1:ceil(opt.N/batch)
        n1 = (i-1)*batch + 1;
        ne = min(opt.N, i*batch);

        if opt.verbose, before = plotBatch(i, batch, opt.N, 50, before); end;
    
        % Compute subjects grad/hess w.r.t. initial velocity
        % --------------------------------------------------
        parfor (n=n1:ne, double(par))
            dat(n) = oneStepGradHessVelocity(dat(n), model, opt);
        end
        
        for n=n1:ne
            
            % Add individual contributions
            % ----------------------------
            for k=1:opt.K
                model.gw(:,:,:,:,k) = model.gw(:,:,:,:,k) + numeric(dat(n).gv) * single(dat(n).z(k));
                model.hw(:,:,:,:,k) = model.hw(:,:,:,:,k) + numeric(dat(n).hv) * single(dat(n).z(k))^2;
            end
            
            % Clear individual grad/hess
            % --------------------------
            dat(n).gv = rmarray(dat(n).gv);
            dat(n).hv = rmarray(dat(n).hv);
        end
        
    end
    if opt.verbose, plotBatchEnd; end;

end

% -----------------
% Residual variance
% -----------------

function [dat, model] = batchGradHessSigma(dat, model, opt)

    % --- Detect parallelisation scheme
    if strcmpi(opt.loop, 'subject') && opt.par > 0
        par = opt.par;
        batch = opt.batch;
    else
        par   = 0;
        batch = 1;
    end

    % Init gradient/hessian
    % ---------------------
    model.gs = 0;
    model.hs = 0;
    
    % --- Batch processing
    if opt.verbose, before = plotBatchBegin('GH Sigma'); end;
    for i=1:ceil(opt.N/batch)
        n1 = (i-1)*batch + 1;
        ne = min(opt.N, i*batch);

        if opt.verbose, before = plotBatch(i, batch, opt.N, 50, before); end;
    
        % Compute subjects grad/hess w.r.t. initial velocity
        % --------------------------------------------------
        parfor (n=n1:ne, double(par))
            dat(n) = oneStepGradHessVelocity(dat(n), model, opt);
        end
        
        for n=n1:ne
            
            % Add individual contributions
            % ----------------------------
            model.gs = model.gs + dat(n).gv(:)' * dat(n).r(:);
            hr = pointwise3(numeric(dat(n).hv), numeric(dat(n).r));
            model.hs = model.hs + dat(n).r(:)' * hr(:);
            clear hr
            
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
    if strcmpi(opt.loop, 'subject')
        loop = '';
        par  = 0;
    else
        loop = opt.loop;
        par  = opt.par;
    end
    
    % --- Transform
    if isfield(todo, 'v')
        if isfield(dat, 'z') && isfield(dat, 'r')
            dat.v = reconstructVelocity('latent', dat.z, 'subspace', model.w, ...
                'residual', dat.r, 'sigma', model.sigma, ...
                'debug', opt.debug, 'output', dat.v, ...
                'loop', loop, 'par', par);
        elseif isfield(dat, 'z')
            dat.v = reconstructVelocity('latent', dat.z, 'subspace', model.w, ...
                'debug', opt.debug, 'output', dat.v, ...
                'loop', loop, 'par', par);
        elseif isfield(dat, 'r') && isfield(model, 'sigma')
            dat.v = reconstructVelocity('residual', dat.r, 'sigma', model.sigma, ...
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
        else
            error('Only affine is not implemented yet')
        end
    end
    if isfield(toclean, 'iphi')
        dat.iphi = rmarray(dat.iphi);
    end
    
    % --- Push
    if isfield(todo, 'pf') || isfield(todo, 'c')
        [dat.pf, dat.c] = pushImage(dat.ipsi, dat.f, opt.lat, ...
            'loop', loop, 'par', par, ...
            'output', {dat.pf, dat.c}, 'debug', opt.debug);
    end
    if isfield(toclean, 'ipsi')
        dat.ipsi = rmarray(dat.ipsi);
    end
    if isfield(todo, 'llm')
        dat.llm = llMatching(opt.model, model.mu, dat.pf, dat.c, ...
            'loop', loop, 'par', par, ...
            'debug', opt.debug);
    end
    
    % --- G/H affine
    if isfield(todo, 'gq') && isfield(todo, 'hq')
        [dat.gq, dat.hq] = ghMatchingAffine(opt.model, ...
            model.mu, dat.pf, dat.c, model.gmu, dat.A, ...
            opt.affine_basis, dat.phi, dat.jac, ...
            'loop', loop, 'par', par);

        [gq, hq] = ghPriorAffine(dat.q, model.regq);
        dat.gq = dat.gq + gq;
        dat.hq = dat.hq + hq;
        clear gq hq
    elseif isfield(todo, 'hq')
        [~, dat.hq] = ghMatchingAffine(opt.model, ...
            model.mu, dat.pf, dat.c, model.gmu, dat.A, ...
            opt.affine_basis, dat.phi, dat.jac, ...
            'loop', loop, 'par', par);

        [~, hq] = ghPriorAffine(dat.q, model.regq);
        dat.hq = dat.hq + hq;
        clear hq
    elseif isfield(todo, 'gq')
        dat.gq = ghMatchingAffine(opt.model, ...
            model.mu, dat.pf, dat.c, model.gmu, dat.A, ...
            opt.affine_basis, dat.phi, dat.jac, ...
            'loop', loop, 'par', par);

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
        [dat.gz, dat.hz] = ghMatchingLatent(opt.model, ...
            model.mu, dat.pf, dat.c, model.gmu, model.w, ...
            'loop', loop, 'par', par);

        [gz, hz] = ghPriorLatent(dat.z, model.regz);
        dat.gz = dat.gz + gz;
        dat.hz = dat.hz + hz;
        clear gz hz
    elseif isfield(todo, 'hz')
        [~, dat.hz] = ghMatchingLatent(opt.model, ...
            model.mu, dat.pf, dat.c, model.gmu, model.w, ...
            'loop', loop, 'par', par);

        [~, hz] = ghPriorLatent(dat.z, model.regz);
        dat.hz = dat.hz + hz;
        clear hz
    elseif isfield(todo, 'gz')
        dat.gz = ghMatchingLatent(opt.model, ...
            model.mu, dat.pf, dat.c, model.gmu, model.w, ...
            'loop', loop, 'par', par);

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
    % ?
    
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
        par = opt.par;
        batch = opt.batch;
    else
        par   = 0;
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
    
    % --- Batch processing
    if opt.verbose, before = plotBatchBegin('Update'); end;
    for i=1:ceil(numel(dat)/batch)
        n1 = (i-1)*batch + 1;
        ne = min(numel(dat), i*batch);
        
        if opt.verbose, before = plotBatch(i, batch, numel(dat), 50, before); end;
        
        parfor (n=n1:ne, double(par))
            dat(n) = oneUpdate(dat(n), model, opt, stodo, clean);
        end
    end
    if opt.verbose, plotBatchEnd; end;
    
end