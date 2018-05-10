function result = lsAffine(model, dq, q0, llm0, mu, f, varargin)
%__________________________________________________________________________
%
% Performs a line search along a direction to find better affine
% parameters. The line search direction is usually found by Gauss-Newton
%
%--------------------------------------------------------------------------
%
% FORMAT result = lsAffine(model, dq, q0, llm0, mu, f, ...)
%
% REQUIRED
% --------
% model - Structure with fields:
%           * 'name'    : 'normal', 'laplace', 'bernoulli' or 'categorical'
%           * ('sigma2'): Normal variance  [1]
%           * ('b')     : Laplace variance [1]
% dq    - Line search direction (ascent if maximising, descent if minimising)
% q0    - Previous parameter value
% llm0  - Previous log-likelihood (matching term)
% mu    - Template (in native space)
% f     - Image (in native space)
%
% KEYWORD ARGUMENTS
% -----------------
% llq0  - Previous log-likelihood (prior term) [compute]
% B     - Affine basis [affine_basis('affine')]
% regq  - Precision matrix of the affine parameters [none]
% rind  - Indices of regularised affine parameters []
% iphi  - Diffeomorphic part [identity]
% Mf    - Image voxel-to-world mappinf [eye(4)]
% Mmu   - Template voxel-ti-world mapping [eye(4)]
% nit   - Number of line-search iterations [6]
% loop  - How to split processing [auto]
% par   - If true, parallelise processing [false]
% itrp  - Template interpolation order [1]
% bnd   - Template Boundary conditions (0=circulant)/[1=neumann]
% 
% OUTPUT
% ------
% Structure with fields:
%   ok   - True if a better parameter value was found
% And if ok == true:
%   q     - New parameter value
%   llm   - New log-likelihood (matching term)
%   llq   - New log-likelihood (prior term)
%   A     - New affine transform
%   pf    - New pushed image
%   c     - New pushed voxel count
%   bb    - Bounding box of pushed image in template space.
%   ipsi  - New complete affine+diffeomorphic mapping
%__________________________________________________________________________

    % --- Parse inputs
    p = inputParser;
    p.FunctionName = 'lsAffine';
    p.addRequired('model',  @(X) isstruct(X) && isfield(X, 'name'));
    p.addRequired('dq',     @checkarray);
    p.addRequired('q0',     @checkarray);
    p.addRequired('llm0',   @isscalar);
    p.addRequired('mu',     @checkarray);
    p.addRequired('f',      @checkarray);
    p.addParameter('llq0',     nan,    @iscalar);
    p.addParameter('B',        [],     @checkarray);
    p.addParameter('rind',     [],     @isnumeric);
    p.addParameter('regq',     []);
    p.addParameter('iphi',     [],     @checkarray);
    p.addParameter('Mf',       eye(4), @(X) isnumeric(X) && issame(size(X), [4 4]));
    p.addParameter('Mmu',      eye(4), @(X) isnumeric(X) && issame(size(X), [4 4]));
    p.addParameter('nit',      6,      @isscalar);
    p.addParameter('match', 'pull', @ischar);
    p.addParameter('par',      false,  @isscalar);
    p.addParameter('loop',     '',     @(X) ischar(X) && any(strcmpi(X, {'slice', 'component', 'none', ''})));
    p.addParameter('itrp',     1,      @isnumeric);
    p.addParameter('bnd',      1,      @(X) isscalar(X) && isnumeric(X));
    p.addParameter('pf',  nan, @(X) isnumeric(X) || isa(X, 'file_array'));
    p.addParameter('c',   nan, @(X) isnumeric(X) || isa(X, 'file_array'));
    p.addParameter('wa',  nan, @(X) isnumeric(X) || isa(X, 'file_array'));
    p.addParameter('wmu', nan, @(X) isnumeric(X) || isa(X, 'file_array'));
    p.addParameter('verbose',  false,  @isscalar);
    p.addParameter('debug',    false,  @isscalar);
    p.parse(model, dq, q0, llm0, mu, f, varargin{:});
    llq0    = p.Results.llq0;
    B       = p.Results.B;
    rind    = p.Results.rind;
    regq    = p.Results.regq;
    iphi    = p.Results.iphi;
    Mf      = p.Results.Mf;
    Mmu     = p.Results.Mmu;
    nit     = p.Results.nit;
    par     = p.Results.par;
    loop    = p.Results.loop;
    verbose = p.Results.verbose;
    debug   = p.Results.debug;
    itrp    = p.Results.itrp;
    bnd     = p.Results.bnd;
    
    matchmode = p.Results.match;
    pf      = p.Results.pf;
    c       = p.Results.c;
    wa      = p.Results.wa;
    wmu     = p.Results.wmu;
    
    if debug, fprintf('* lsAffine\n'); end
    
    cat = any(strcmpi(model.name, {'bernoulli', 'binomial', 'categorical', 'multinomial'}));
    
    % --- Save previous pushed/pull
    pf_bck = false;
    if isa(pf, 'file_array') && exist(pf.fname, 'file')
        pf_bck = true;
        [path, name, ext] = fileparts(pf.fname);
        copyfile(pf.fname, fullfile(path, [name '_bck' ext]));
    end
    c_bck = false;
    if isa(c, 'file_array') && exist(c.fname, 'file')
        c_bck = true;
        [path, name, ext] = fileparts(c.fname);
        copyfile(c.fname, fullfile(path, [name '_bck' ext]));
    end
    wa_bck = false;
    if isa(wa, 'file_array') && exist(wa.fname, 'file')
        wa_bck = true;
        [path, name, ext] = fileparts(wa.fname);
        copyfile(wa.fname, fullfile(path, [name '_bck' ext]));
    end
    wmu_bck = false;
    if isa(wmu, 'file_array') && exist(wmu.fname, 'file')
        wmu_bck = true;
        [path, name, ext] = fileparts(wmu.fname);
        copyfile(wmu.fname, fullfile(path, [name '_bck' ext]));
    end
    
    % --- Set some default parameter value
    if isempty(B)
        error('Missing affine basis')
    end
    if ~isempty(regq) && isempty(rind)
        rind = (1:size(regq,1)) + numel(q0) - size(regq, 1);
    end
    if isnan(llq0)
        if ~isempty(regq),  llq0 = llPriorAffine(q0(rind), regq, 'fast', 'debug', debug);
        else,               llq0 = 0; end
    end
    dim  = [size(mu) 1 1];
    lat = dim(1:3);
    if isempty(iphi)
        iphi = spm_warps('identity', lat);
    end
    
    % --- Load some data (in case it is on disk)
    q0 = numeric(q0);
    dq = numeric(dq);
    
    % --- Initialise line search
    armijo = 1;           % Armijo facto
    ll0    = llm0 + llq0; % Log-likelihood (only parts that depends on q)
    dimf   = [size(f) 1 1];
    latf   = dimf(1:3);
    dimmu   = [size(mu) 1 1];
    latmu   = dimmu(1:3);
    
    if verbose
        printInfo('header');
        printInfo('initial', ll0, llm0, llq0);
    end
    
    % --- Loop
    for i=1:nit
        q = q0 + dq / armijo;
        A = exponentiateAffine(q, B, 'debug', debug);
        ipsi = reconstructIPsi(A, iphi, 'lat', latf, 'Mf', Mf, 'Mmu', Mmu, 'debug', debug);
        if strcmpi(matchmode, 'push')
            [pf, c, bb] = pushImage(ipsi, f, latmu, 'circ', ~bnd, 'par', par, 'loop', loop, 'debug', debug, 'output', {pf, c});
            llm = llMatching(model, mu, pf, c, 'bb', bb, 'par', par, 'loop', loop, 'debug', debug);
        elseif strcmpi(matchmode, 'pull')
            if cat
                wa = pullTemplate(ipsi, mu, 'order', itrp, 'bnd', ~bnd, 'par', par, 'output', wa, 'debug', debug);
                wmu = reconstructProbaTemplate(wa, 'output', wmu, 'loop', loop, 'par', par, 'debug', debug);
                if ~isa(wa, 'file_array')
                    clear wa
                end
            else
                wmu = pullTemplate(ipsi, mu, 'order', itrp, 'bnd', ~bnd, 'par', par, 'output', wmu, 'debug', debug);
            end
            llm = llMatching(model, wmu, f, 'par', par, 'loop', loop, 'debug', debug);
        end
        if checkarray(regq)
            llq = llPriorAffine(q(rind), regq, 'fast', 'debug', debug);
        else
            llq = 0;
        end
        ll = llm + llq;
        
        if verbose, printInfo(i, ll0, llm, llq); end
        
        if ll <= ll0
            if verbose, printInfo('failed'); end
            armijo = armijo * 2;
        else
            if verbose, printInfo('success'); end
            result     = struct;
            result.ok  = true;
            result.q   = q;
            result.A   = A;
            if checkarray(regq)
                result.llq = llPriorAffine(q(rind), regq, 'debug', debug);
            end
            result.llm = llm;
            result.ipsi = ipsi;
            if strcmpi(matchmode, 'pull')
                [pf, c, bb] = pushImage(ipsi, f, latmu, 'circ', ~bnd, 'par', par, 'loop', loop, 'debug', debug, 'output', {pf, c});
                result.wmu = wmu;
                result.wa  = wa;
            end
            result.pf = pf;
            result.c  = c;
            result.bb = bb;
            return
        end
    end
    
    if verbose, printInfo('end'); end
    
    result       = struct;
    result.ok    = false;
    if pf_bck
        [path, name, ext] = fileparts(pf.fname);
        movefile(fullfile(path, [name '_bck' ext]), pf.fname);
        
    end
    if c_bck
        [path, name, ext] = fileparts(c.fname);
        movefile(fullfile(path, [name '_bck' ext]), c.fname);
    end
    if wa_bck
        [path, name, ext] = fileparts(wa.fname);
        movefile(fullfile(path, [name '_bck' ext]), wa.fname);
    end
    if wmu_bck
        [path, name, ext] = fileparts(wmu.fname);
        movefile(fullfile(path, [name '_bck' ext]), wmu.fname);
    end

end

function printInfo(which, oll, llm, llq)
    if ischar(which) 
        if strcmpi(which, 'header')
            fprintf('Q - LineSearch | Armijo  | %12s = %12s + %12s | %12s\n', 'RLL', 'LL-Match', 'RLL-Prior', 'LL-Diff');
        elseif strcmpi(which, 'initial')
            fprintf('Q - LineSearch | Initial | %12.6f = %12.6f + %12.6f \n', oll, llm, llq);
        elseif strcmpi(which, 'failed')
            fprintf('| Failed\n');
        elseif strcmpi(which, 'success')
            fprintf('| Success\n');
        elseif strcmpi(which, 'end')
            fprintf('Q - LineSearch | Complete failure\n');
        end
    else
        fprintf('Q - LineSearch | Try %3d | %12.6f = %12.6f + %12.6f | %12.6f ', which, llm+llq, llm, llq, llm+llq-oll);
    end
end
        