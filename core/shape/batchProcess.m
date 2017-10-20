function varargout = batchProcess(id, varargin)
% Collection of functions that require updating some subject-specific
% variables. Depending on the needs, as much as possible is done in order
% to minimise i/o, memory and computational costs.
%
% FORMAT [dat, model] = batchProcess('Init', dat, model, opt)
% > Initialise latent coordinates and insure they are zero-centered (but
%   orthogonalisation is not done here).
%   Update>> dat.z, dat.zz, dat.Sz, model.z, model.zz, model.Sz
%
% FORMAT dat = batchProcess('Rotate', dat, opt, R)
% > Apply a rotation matrix to all latent coordinates
%   Update>> dat.z, dat.zz, dat.Sz
%
% FORMAT dat = batchProcess('Centre', dat, opt, (mean))
% > Subtract the mean to all latent coordinates
%   Update>> dat.z, dat.zz
%
% FORMAT [dat, model] = batchProcess('E', dat, model, opt)
%
% FORMAT [dat, model] = batchProcess('M', dat, model, opt)
%
% FORMAT dat = batchProcess('Update', dat, model, opt, {...}, 'clean', {...})
% > Generic tool to update specified subject variables. Names of variables
%   to update are provided in the first list, and names of variables to
%   discard in the end (to save memory/disk) can be provided in the second
%   list.
%   The model is never updated, consequently, the processing is not really
%   performed by batch but distributed as specified by opt.loop/opt.par.


    switch lower(id)
        case 'init'
            [varargout{1:nargout}] = batchInit(varargin{:});
        case 'rotate'
            [varargout{1:nargout}] = batchRotate(varargin{:});
        case 'centre'
            [varargout{1:nargout}] = batchCentre(varargin{:});
        case 'm0'
            [varargout{1:nargout}] = batchM0(varargin{:});
        case 'em'
            [varargout{1:nargout}] = batchEM(varargin{:});
        case 'e'
            [varargout{1:nargout}] = batchE(varargin{:});
        case 'm'
            [varargout{1:nargout}] = batchM(varargin{:});
        case 'll'
            [varargout{1:nargout}] = batchLL(varargin{:});
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
%    Init
% -------------------------------------------------------------------------

function [dat, model] = batchInit(mode, dat, model, opt)

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
    if opt.verbose, before = plotBatchBegin('Init'); end;
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
    if opt.verbose, before = plotBatchBegin('Center'); end;
    for n=1:opt.N
        if opt.verbose, before = plotBatch(n, 1, opt.N, 50, before); end;
        z = numeric(dat(n).z) - mz/opt.N;
        
        dat(n).z  = saveOnDisk(dat(n).z, z);
        dat(n).zz = saveOnDisk(dat(n).z, z*z');
        dat(n).Sz = saveOnDisk(dat(n).Sz, zeros(opt.K));
    end
    if opt.verbose, plotBatchEnd; end;
    model.z  = saveOnDisk(model.z,  zeros(opt.K, 1));
    model.zz = saveOnDisk(model.zz, mzz - mz*mz'/opt.N);
    model.Sz  = saveOnDisk(model.Sz,  zeros(opt.K));

end

% -------------------------------------------------------------------------
%    Center
% -------------------------------------------------------------------------

function dat = batchCentre(dat, opt, mean)

    % --- Don't parallelise here, there is no need
    
    % Compute the mean if needed
    % --------------------------
    if nargin < 3
        if opt.verbose, before = plotBatchBegin('Accumulate'); end;
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
    if opt.verbose, before = plotBatchBegin('Centre'); end;
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

function dat = batchRotate(dat, opt, R)

    % --- Don't parallelise here, there is no need
    
    % Rotate subjects
    % ---------------
    if opt.verbose, before = plotBatchBegin('Rotate'); end;
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

function dat = oneStepE(dat, model, opt)

    % Detect parallelisation scheme
    % ------------------------
    if strcmpi(opt.loop, 'subject')
        opt.loop = '';
        opt.par  = 0;
    else
        opt.loop = opt.loop;
        opt.par  = opt.par;
    end
    
    % Gauss-Newton iterations
    % -----------------------
    % It is useful to actually find a mode of the posterior (and not only
    % an improved value) when we use the Laplace precision for the update
    % of W. In that case, setting gnit > 1 might help converge faster.
    for i=1:opt.gnit

        % Compute gradient/hessian
        % ------------------------

        [dat.g, dat.h] = ghMatchingLatent(opt.model, ...
            model.mu, dat.pf, dat.c, model.gmu, model.w, ...
            'loop', opt.loop, 'par', opt.par, 'debug', opt.debug, ...
            'output', {dat.g, dat.h});

        [g, h] = ghPriorLatent(dat.z, model.regz, 'debug', opt.debug);
        dat.g = dat.g + g;
        dat.h = dat.h + h;

        dat.h = loadDiag(dat.h); % Additional regularisation for robustness
        
        % Compute search direction
        % ------------------------
        dat.dz = -dat.h \ dat.g;

        % Line search
        % -----------
        [okz, z, llm, ~, v, iphi, pf, c, ipsi] = lsLatent(...
            opt.model, dat.dz, dat.z, dat.v, dat.llm, ...
            model.w, model.mu, dat.f, ...
            'regz', model.regz, ...
            'Mf', dat.Mf, 'Mmu', model.Mmu, 'nit', opt.lsit, ...
            'itgr', opt.itgr, 'prm', opt.prm);

        % Store better values
        % -------------------
        dat.okz = okz;
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
        [~, dat.h] = ghMatchingLatent(opt.model, ...
            model.mu, dat.pf, dat.c, model.gmu, model.w);

        [~, h] = ghPriorLatent(dat.z, model.regz);
        dat.h = dat.h + h;

        dat.h = loadDiag(dat.h); % Additional regularisation for robustness
        
        dat.Sz = inv(dat.h);
    end
        
    % Cleaning
    % --------
    % I should probably clear variables and remove files that are not
    % useful anymore. This will cause less disk and broadband usage.
    toclean = {'g', 'h', 'iphi', 'ipsi'};
    for i=1:numel(toclean)
        field = toclean{i};
        dat.(field) = rmarray(dat.(field));
    end
    if isfield(dat, 'dz')
        dat = rmfield(dat, 'dz');
    end
    if isfield(dat, 'll')
        dat = rmfield(dat, 'll');
    end
end

function [dat, model] = batchE(dat, model, opt)

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
    mz         = zeros(opt.K, 1);
    mzz        = zeros(opt.K, opt.K);
    mS         = zeros(opt.K, opt.K);
    okz        = 0;
    
    % --- Batch processing
    if opt.verbose, before = plotBatchBegin('E'); end;
    for i=1:ceil(opt.N/batch)
        n1 = (i-1)*batch + 1;
        ne = min(opt.N, i*batch);

        if opt.verbose, before = plotBatch(i, batch, opt.N, 50, before); end;
    
        % Compute subjects grad/hess w.r.t. initial velocity
        % --------------------------------------------------
        parfor (n=n1:ne, double(par))
            dat(n) = oneStepE(dat(n), model, opt);
        end
        
        for n=n1:ne
            
            % Add individual contributions
            % ----------------------------
            mz         = mz  + dat(n).z;
            mzz        = mzz + dat(n).zz;
            mS         = mS  + dat(n).Sz;
            model.llm  = model.llm  + dat(n).llm;
            okz        = okz + dat(n).okz;
            
        end
        
    end
    if opt.verbose, fprintf(' | %4d / %4d', okz, opt.N); end;
    if opt.verbose, plotBatchEnd; end;
    
    model.z    = saveOnDisk(model.z, mz);
    model.zz   = saveOnDisk(model.zz, mzz);
    model.Sz    = saveOnDisk(model.Sz, mS);

end

% -------------------------------------------------------------------------
%    M-step
% -------------------------------------------------------------------------

function dat = oneStepM(dat, model, opt)

    [dat.gv, dat.hv] = ghMatchingVel(opt.model, ...
        model.mu, dat.pf, dat.c, model.gmu, ...
        'output', {dat.gv, dat.hv});
end

function [dat, model] = batchM(dat, model, opt)

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
    if opt.verbose, before = plotBatchBegin('M'); end;
    for i=1:ceil(opt.N/batch)
        n1 = (i-1)*batch + 1;
        ne = min(opt.N, i*batch);

        if opt.verbose, before = plotBatch(i, batch, opt.N, 50, before); end;
    
        % Compute subjects grad/hess w.r.t. initial velocity
        % --------------------------------------------------
        parfor (n=n1:ne, double(par))
            dat(n) = oneStepM(dat(n), model, opt);
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
    
    if isfield(todo, 'v')
        dat.v = reconstructVelocity('latent', dat.z, 'subspace', model.w, ...
            'debug', opt.debug, 'output', dat.v, ...
            'loop', loop, 'par', par);
    end
    if isfield(todo, 'iphi')
        dat.iphi = exponentiateVelocity(dat.v, 'iphi', ...
            'itgr', opt.itgr, 'vs', opt.vs, ...
            'prm', opt.prm, 'debug', opt.debug, 'output', dat.iphi);
    end
    if isfield(toclean, 'v')
        dat.v = rmarray(dat.v);
    end
    if isfield(todo, 'ipsi')
        latf = [size(dat.f) 1];
        latf = latf(1:3);
        dat.ipsi = reconstructIPsi(eye(4), dat.iphi, ...
            'lat', latf, 'Mf', dat.Mf, 'Mmu', model.Mmu, ...
            'output', dat.ipsi, 'debug', opt.debug);
    end
    if isfield(toclean, 'iphi')
        dat.iphi = rmarray(dat.iphi);
    end
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
    if isfield(todo, 'g') && isfield(todo, 'h')
        [dat.g, dat.h] = ghMatchingLatent(opt.model, ...
            model.mu, dat.pf, dat.c, model.gmu, model.w, ...
            'loop', loop, 'par', par);

        [g, h] = ghPriorLatent(dat.z, model.regz);
        dat.g = dat.g + g;
        dat.h = dat.h + h;
    elseif isfield(todo, 'h')
        [~, dat.h] = ghMatchingLatent(opt.model, ...
            model.mu, dat.pf, dat.c, model.gmu, model.w, ...
            'loop', loop, 'par', par);

        [~, h] = ghPriorLatent(dat.z, model.regz);
        dat.h = dat.h + h;
    elseif isfield(todo, 'g')
        dat.g = ghMatchingLatent(opt.model, ...
            model.mu, dat.pf, dat.c, model.gmu, model.w, ...
            'loop', loop, 'par', par);

        g = ghPriorLatent(dat.z, model.regz);
        dat.g = dat.g + g;
    end
    if isfield(todo, 'h')
        dat.h = loadDiag(dat.h);
    end
    
    if isfield(toclean, 'pf')
        dat.pf = rmarray(dat.pf);
    end
    if isfield(toclean, 'c')
        dat.c = rmarray(dat.c);
    end
    if isfield(todo, 'lllz')
        dat.lllz = llLaplace(dat.h);
    end
    if isfield(todo, 'S')
        dat.Sz = inv(dat.h);
    end
    if isfield(toclean, 'g')
        dat.g = rmarray(dat.g);
    end
    if isfield(toclean, 'h')
        dat.h = rmarray(dat.h);
    end
    if isfield(todo, 'llz')
        dat.llz = llPriorLatent(dat.z, model.regz, 'debug', opt.debug);
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
    for i=1:ceil(opt.N/batch)
        n1 = (i-1)*batch + 1;
        ne = min(opt.N, i*batch);
        
        if opt.verbose, before = plotBatch(i, batch, opt.N, 50, before); end;
        
        parfor (n=n1:ne, double(par))
            dat(n) = oneUpdate(dat(n), model, opt, stodo, clean);
        end
    end
    if opt.verbose, plotBatchEnd; end;
    
end