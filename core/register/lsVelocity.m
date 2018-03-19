function result = lsVelocity(model, dv, r0, match0, mu, f, varargin)
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
% bnd  - Differential operator boundary conditions: [0]/1/2/3$
% geod - eight of the geodesic prior [0]
% reg0 - Previous prior term [compute]
% vel0 - Previous geodesic term [compute]
%
% A     - Affine transform [eye(4)]
% Mf    - Image voxel-to-world mappinf [eye(4)]
% Mmu   - Template voxel-to-world mapping [eye(4)]
%
% match  - Mode for computing thr matching term: ['pull']/'push'
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
% And if ok == true:
%     v     - New initial velocity
%     match - New log-likelihood (matching term)
%     pf    - New pushed image 
%     c     - New pushed voxel count
%     wmu   - New pulled template
%    ipsi   - New complete affine+diffeomorphic mapping
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
    p.addParameter('vs',   [],      @isnumeric);
    p.addParameter('itgr', nan,     @isscalar);
    p.addParameter('bnd',  0,       @(X) isscalar(X) && isnumeric(X));
    p.addParameter('geod', 0,       @isscalar);
    p.addParameter('reg0', nan,     @isscalar);
    p.addParameter('vel0', nan,     @isscalar);
    
    p.addParameter('A',     eye(4), @(X) isnumeric(X) && issame(size(X), [4 4]));
    p.addParameter('Mf',    eye(4), @(X) isnumeric(X) && issame(size(X), [4 4]));
    p.addParameter('Mmu',   eye(4), @(X) isnumeric(X) && issame(size(X), [4 4]));
    
    p.addParameter('match', 'pull', @ischar);
    
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
    geod    = p.Results.geod;
    reg0    = p.Results.reg0;
    vel0    = p.Results.vel0;
    
    A       = p.Results.A;
    Mf      = p.Results.Mf;
    Mmu     = p.Results.Mmu;
    
    nit     = p.Results.nit;
    par     = p.Results.par;
    loop    = p.Results.loop;
    
    matchmode = p.Results.match;
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
        reg0 = (1-geod) * llPriorVelocity(r0, 'fast', 'vs', vs, 'prm', prm, 'bnd', bnd, 'debug', debug);

    end 
    if isnan(vel0)
        if geod
            vel0 = geod * llPriorVelocity(v0, 'fast', 'vs', vs, 'prm', prm, 'bnd', bnd, 'debug', debug);
        else
            vel0 = 0;
        end
    end
    
    % --- Initialise line search
    armijo = 1;             % Armijo factor
    ll0    = match0 + reg0 + vel0; % Log-likelihood (only parts that depends on v)
    dimf   = [size(f) 1 1];
    latf   = dimf(1:3);
    dimmu  = [size(mu) 1 1];
    latmu  = dimmu(1:3);
    
    if verbose
        printInfo('header');
        printInfo('initial', ll0, match0, reg0, vel0);
    end
    
    % --- Loop
    for i=1:nit
        v    = single(v0 + dv / armijo);
        iphi = exponentiateVelocity(v, 'iphi', 'itgr', itgr, 'vs', vs, 'prm', prm, 'bnd', bnd, 'debug', debug);
        ipsi = reconstructIPsi(A, iphi, 'lat', latf, 'Mf', Mf, 'Mmu', Mmu, 'debug', debug);
        if strcmpi(matchmode, 'push')
            [pf, c, bb] = pushImage(ipsi, f, latmu, 'par', par, 'loop', loop, 'debug', debug, 'output', {pf, c});
            match = llMatching(model, mu, pf, c, 'bb', bb, 'par', par, 'loop', loop, 'debug', debug);
        elseif strcmpi(matchmode, 'pull')
            if cat
                wa = pullTemplate(ipsi, mu, 'par', par, 'output', wa, 'debug', debug);
                wmu = reconstructProbaTemplate(wa, 'output', wmu, 'loop', loop, 'par', par, 'debug', debug);
            else
                wmu = pullTemplate(ipsi, mu, 'par', par, 'output', wmu, 'debug', debug);
            end
            match = llMatching(model, wmu, f, 'par', par, 'loop', loop, 'debug', debug);
        end
        r = r0 + dv / armijo;
        reg = (1-geod) * llPriorVelocity(r,  'fast', 'vs', vs, 'prm', lam*prm, 'bnd', bnd, 'debug', debug);
        if geod
            vel = geod * llPriorVelocity(v, 'fast', 'vs', vs, 'prm', prm, 'bnd', bnd, 'debug', debug);
        else
            vel = 0;
        end
        ll  = match + reg + vel;
        
        if verbose, printInfo(i, ll0, match, reg, vel); end
        
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
            result.iphi = iphi;
            result.ipsi = ipsi;
            if strcmpi(matchmode, 'pull')
                [pf, c, bb] = pushImage(ipsi, f, latmu, 'par', par, 'loop', loop, 'debug', debug, 'output', {pf, c});
                result.wmu = wmu;
                result.wa = wa;
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

function printInfo(which, oll, llm, llr, llv)
    if ischar(which) 
        if strcmpi(which, 'header')
            fprintf('R - LineSearch | Armijo  | %12s = %12s + %12s | %12s\n', 'RLL', 'LL-Match', 'RLL-Prior', 'LL-Diff');
        elseif strcmpi(which, 'initial')
            fprintf('R - LineSearch | Initial | %12.4f = %12.4f + %12.4f + %12.4f \n', oll, llm, llr, llv);
        elseif strcmpi(which, 'failed')
            fprintf('| Failed\n');
        elseif strcmpi(which, 'success')
            fprintf('| Success\n');
        elseif strcmpi(which, 'end')
            fprintf('R - LineSearch | Complete failure\n');
        end
    else
        fprintf('R - LineSearch | Try %3d | %12.4f = %12.4f + %12.4f + %12.4f | %12.4f ', which, llm+llr+llv, llm, llr, llv, llm+llr+llv-oll);
    end
end
        