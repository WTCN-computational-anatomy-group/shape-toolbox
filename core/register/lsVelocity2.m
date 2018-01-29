function result = lsVelocity2(model, dv, r0, match0, mu, f, varargin)
%__________________________________________________________________________
% 
% Performs a line search along a direction to find better initial
% velocity. The line search direction is usually found by Gauss-Newton.
%
%--------------------------------------------------------------------------
%
% FORMAT [ok, ...] = lsVelocity(model, dv, r0, match0, mu, f, ...)
%
% REQUIRED
% --------
% model  - Structure with fields:
%           * 'name'    : 'normal', 'laplace', 'bernoulli' or 'categorical'
%           * ('sigma2'): Normal variance  [1]
%           * ('b')     : Laplace variance [1]
% dv     - Line search direction (ascent if maximising, descent if minimising)
% r0     - Previous residual field
% match0 - Previous log-likelihood (matching term)
% mu     - Template (in native space)
% f      - Image (in native space)
%
% KEYWORD ARGUMENTS
% -----------------
% v0   - Previous velocity field [r0]
% lam  - Precision modulation [1]
% prm  - Regularisation parameters (L matrix) [0.0001 0.001 0.2 0.05 0.2]
% vs   - Lattice voxel size [[] = from Mmu]
% itgr - Number of integration steps for geodesic shooting [NaN = auto]
% bnd  - Differential operator boundary conditions: [0]/1/2/3
%
% reg0  - Previous prior term [NaN = compute]
% A     - Affine transform [eye(4)]
% Mf    - Image voxel-to-world mappinf [eye(4)]
% Mmu   - Template voxel-to-world mapping [eye(4)]
%
% match  - Mode for computing thr matching term: ['pull']/'push'
% itrp   - Template interpolation order (if match is 'pull') [1]
% tplbnd - Template boundary condition (if match is 'pull') [1]
%
% nit  - Number of line-search iterations [6]
% loop - How to split processing [auto]
% par  - If true, parallelise processing [false]
%
% pf   - file array where to store output pushed image
% c    - file array where to store output pushed count
% wa   - file array where to store output pulled log-template
% wmu  - file array where to store output pulled template
%
% OUTPUT
% ------
% Structure with fields:
%     ok    - True if a better parameter value was found
%     v     - New initial velocity
%     match - New log-likelihood (matching term)
%     pf    - New pushed image 
%     c     - New pushed voxel count
%     wmu   - New pulled template
%
%--------------------------------------------------------------------------
% 
% The objective function is supposed to be of the form
%   E = logMatch(f|mu,v) + logN(v | mean, inv(lam*L))
%
% The regularisation is split in two parts (lam and L), because L (alone)
% is used to exponentiate the velocity field (i.e., it is also
% parameterises shooting.)
%
% The classical diffeomorphic fitting is done by setting 
%    mean = NaN and lam = 1
%__________________________________________________________________________

    % ---------------------------------------------------------------------
    % Parse inputs
    % ------------
    dft_prm = [0.0001 0.001 0.2 0.05 0.2];
    
    p = inputParser;
    p.FunctionName = 'lsVelocity';
    p.addRequired('model',  @(X) isstruct(X) && isfield(X, 'name'));
    p.addRequired('dv',     @checkarray);
    p.addRequired('r0',     @checkarray);
    p.addRequired('match0', @isscalar);
    p.addRequired('mu',     @checkarray);
    p.addRequired('f',      @checkarray);
    
    p.addParameter('v0',   [],      @checkarray);
    p.addParameter('lam',  1,       @isscalar);
    p.addParameter('prm',  dft_prm, @(X) length(X) == 5);
    p.addParameter('pgprior', false,@(X) isscalar(X) && islogical(X));
    p.addParameter('vs',   [],      @isnumeric);
    p.addParameter('itgr', nan,     @isscalar);
    p.addParameter('bnd',  0,       @(X) isscalar(X) && isnumeric(X));
    
    p.addParameter('reg0',  nan,    @isscalar);
    p.addParameter('A',     eye(4), @(X) isnumeric(X) && issame(size(X), [4 4]));
    p.addParameter('Mf',    eye(4), @(X) isnumeric(X) && issame(size(X), [4 4]));
    p.addParameter('Mmu',   eye(4), @(X) isnumeric(X) && issame(size(X), [4 4]));
    
    p.addParameter('match', 'pull', @ischar);
    p.addParameter('itrp', 1, @(X) isscalar(X) && isnumeric(X) );
    p.addParameter('tplbnd', 1, @(X) isscalar(X) && isnumeric(X) );
    
    p.addParameter('nit',  6,     @isscalar);
    p.addParameter('par',  false, @isscalar);
    p.addParameter('loop', '',    @(X) ischar(X) && any(strcmpi(X, {'slice', 'component', 'none', ''})));
    
    p.addParameter('pf',  nan, @(X) isnumeric(X) || isa(X, 'file_array'));
    p.addParameter('c',   nan, @(X) isnumeric(X) || isa(X, 'file_array'));
    p.addParameter('wa',  nan, @(X) isnumeric(X) || isa(X, 'file_array'));
    p.addParameter('wmu', nan, @(X) isnumeric(X) || isa(X, 'file_array'));
    
    p.addParameter('verbose', false, @isscalar);
    p.addParameter('debug',   false, @isscalar);
    
    p.parse(model, dv, r0, match0, mu, f, varargin{:});
    
    v0      = p.Results.v0;
    lam     = p.Results.lam;
    prm     = p.Results.prm;
    vs      = p.Results.vs;
    itgr    = p.Results.itgr;
    bnd     = p.Results.bnd;
    pgprior = p.Results.pgprior;
    
    reg0    = p.Results.reg0;
    A       = p.Results.A;
    Mf      = p.Results.Mf;
    Mmu     = p.Results.Mmu;
    
    matchmode = p.Results.match;
    itrp      = p.Results.itrp;
    tplbnd    = p.Results.tplbnd;
    
    nit     = p.Results.nit;
    par     = p.Results.par;
    loop    = p.Results.loop;
    
    pf      = p.Results.pf;
    c       = p.Results.c;
    wa      = p.Results.wa;
    wmu     = p.Results.wmu;
    
    verbose = p.Results.verbose;
    debug   = p.Results.debug;
    % ---------------------------------------------------------------------
    
    if debug, fprintf('* lsVelocity\n'); end
    
    cat = any(strcmpi(model.name, {'bernoulli', 'binomial', 'categorical', 'multinomial'}));
    
    % --- Save previous pushed/pull
    if ~isa(pf, 'file_array')
        [path, name, ext] = fileparts(pf.fname);
        copyfile(pf.fname, fullfile(path, [name '_bck' ext]));
    end
    if ~isa(c, 'file_array')
        [path, name, ext] = fileparts(c.fname);
        copyfile(c.fname, fullfile(path, [name '_bck' ext]));
    end
    if ~isa(wa, 'file_array')
        [path, name, ext] = fileparts(wa.fname);
        copyfile(wa.fname, fullfile(path, [name '_bck' ext]));
    end
    if ~isa(wmu, 'file_array')
        [path, name, ext] = fileparts(wmu.fname);
        copyfile(wmu.fname, fullfile(path, [name '_bck' ext]));
    end
    
    % --- Load some data (in case it is on disk)
    r0 = numeric(r0);
    v0 = numeric(v0);
    dv = numeric(dv);
    if ~checkarray(v0)
        v0 = r0;
    end
    
    % --- Set some default parameter value
    if ~checkarray(vs)
        vs = sqrt(sum(Mmu(1:3,1:3).^2)); 
    end
    if isnan(reg0)
        reg0 = llPriorVelocity(r0, 'fast', 'vs', vs, 'prm', lam*prm, 'bnd', bnd, 'debug', debug);
        if pgprior
            reg0 = reg0 + llPriorVelocity(v0, 'fast', 'vs', vs, 'prm', prm, 'bnd', bnd, 'debug', debug);
        end
    end 
    
    % --- Initialise line search
    armijo = 1;             % Armijo factor
    ok     = false;         % Found a better ll ?
    ll0    = match0 + reg0; % Log-likelihood (only parts that depends on v)
    dimf   = [size(f) 1 1];
    latf   = dimf(1:3);
    dimmu  = [size(mu) 1 1];
    latmu  = dimmu(1:3);
    
    if verbose
        printInfo('header');
        printInfo('initial', ll0, match0, reg0);
    end
    
    % --- Loop
    for i=1:nit
        v    = single(v0 + dv / armijo);
        iphi = exponentiateVelocity(v, 'iphi', 'itgr', itgr, 'vs', vs, 'prm', prm, 'bnd', bnd, 'debug', debug);
        ipsi = reconstructIPsi(A, iphi, 'lat', latf, 'Mf', Mf, 'Mmu', Mmu, 'debug', debug);
        clear iphi
        if strcmpi(matchmode, 'push')
            [pf, c, bb] = pushImage(ipsi, f, latmu, 'par', par, 'loop', loop, 'debug', debug, 'output', {pf, c});
            match = llMatching(model, mu, pf, c, 'bb', bb, 'par', par, 'loop', loop, 'debug', debug);
        elseif strcmpi(matchmode, 'pull')
            if cat
                wa = pullTemplate(ipsi, mu, 'par', par, 'output', wa, 'debug', debug);
                wmu = reconstructProbaTemplate(wa, 'output', wmu, 'loop', loop, 'par', par, 'debug', debug);
                if ~isa(wa, 'file_array')
                    clear wa
                end
            else
                wmu = pullTemplate(ipsi, mu, 'par', par, 'output', wmu, 'debug', debug);
            end
            match = llMatching(model, wmu, f, 'par', par, 'loop', loop, 'debug', debug);
        end
        r = r0 + dv / armijo;
        reg = llPriorVelocity(r,  'fast', 'vs', vs, 'prm', lam*prm, 'bnd', bnd, 'debug', debug);
        if pgprior
            reg = reg + llPriorVelocity(v,  'fast', 'vs', vs, 'prm', prm, 'bnd', bnd, 'debug', debug);
        end
        ll  = match + reg;
        
        if verbose, printInfo(i, ll0, match, reg); end
        
        if ll <= ll0
            if verbose, printInfo('failed'); end
            armijo = armijo * 2;
        else
            if verbose, printInfo('success'); end
            ok  = true;
            result = struct;
            result.ok = ok;
            result.match = match;
            result.v = v;
            result.r = r;
            if strcmpi(matchmode, 'pull')
                [pf, c, bb] = pushImage(ipsi, f, latmu, 'par', par, 'loop', loop, 'debug', debug, 'output', {pf, c});
                result.wmu = wmu;
            end
            result.pf = pf;
            result.c  = c;
            result.bb = bb;
            return
        end
    end
    
    if verbose, printInfo('end'); end
    
    result       = struct;
    result.ok    = ok;
    result.match = match0;
    result.v     = v0;
    result.r     = r0;
    if ~isa(pf, 'file_array')
        [path, name, ext] = fileparts(pf.fname);
        copyfile(pf.fname, fullfile(path, [name(1:end-4) ext]));
        result.pf = pf;
    end
    if ~isa(c, 'file_array')
        [path, name, ext] = fileparts(c.fname);
        copyfile(c.fname, fullfile(path, [name(1:end-4) ext]));
        result.c = c;
    end
    if ~isa(wa, 'file_array')
        [path, name, ext] = fileparts(wa.fname);
        copyfile(wa.fname, fullfile(path, [name(1:end-4) ext]));
        result.wa = wa;
    end
    if ~isa(wmu, 'file_array')
        [path, name, ext] = fileparts(wmu.fname);
        copyfile(wmu.fname, fullfile(path, [name(1:end-4) ext]));
        result.wmu = wmu;
    end

end

function printInfo(which, oll, llm, llr)
    if ischar(which) 
        if strcmpi(which, 'header')
            fprintf('R - LineSearch | Armijo  | %12s = %12s + %12s | %12s\n', 'RLL', 'LL-Match', 'RLL-Prior', 'LL-Diff');
        elseif strcmpi(which, 'initial')
            fprintf('R - LineSearch | Initial | %12.4f = %12.4f + %12.4f \n', oll, llm, llr);
        elseif strcmpi(which, 'failed')
            fprintf('| Failed\n');
        elseif strcmpi(which, 'success')
            fprintf('| Success\n');
        elseif strcmpi(which, 'end')
            fprintf('R - LineSearch | Complete failure\n');
        end
    else
        fprintf('R - LineSearch | Try %3d | %12.4f = %12.4f + %12.4f | %12.4f ', which, llm+llr, llm, llr, llm+llr-oll);
    end
end
        