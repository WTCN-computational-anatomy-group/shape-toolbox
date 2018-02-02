function result = lsLatent(model, dz, z0, v0, llm0, W, mu, f, varargin)
%__________________________________________________________________________
%
% Performs a line search along a direction to find better latent
% coordinates. The line search direction is usually found by Gauss-Newton.
%
%--------------------------------------------------------------------------
%
% FORMAT result = lsLatent(model, dz, z0, v0, llm0, W, mu, f, ...)
%
% REQUIRED
% --------
% model - Structure with fields:
%           * 'name'    : 'normal', 'laplace', 'bernoulli' or 'categorical'
%           * ('sigma2'): Normal variance  [1]
%           * ('b')     : Laplace variance [1]
% dz   - Line search direction (ascent if maximising, descent if minimising)
% z0   - Previous parameter value
% v0   - Previous velocity field
% llm0 - Previous log-likelihood (matching term)
% W    - Principal subspace
% mu   - Template (in native space)
% f    - Image (in native space)
%
% KEYWORD ARGUMENTS
% -----------------
% regz - Precision matrix of the latent parameters [none]
% prm  - Differential operator parameters [0.0001 0.001 0.2 0.05 0.2]
% itgr - Number of integration steps for geodesic shooting [auto]
% geod - Weight of the geodesic prior [0]
% llz0 - Previous log-likelihood (prior term) [compute]
% llv0 - Previous log-likelihood (geodesic term) [compute]
%
% A    - Affine transform [eye(4)]
% Mf   - Image voxel-to-world mappinf [eye(4)]
% Mmu  - Template voxel-to-world mapping [eye(4)]
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
%   ok   - True if a better parameter value was found
% And if ok == true:
%   z    - New parameter value
%   llm  - New log-likelihood (matching term)
%   llz  - New log-likelihood (prior term)
%   v    - New velocity field
%   iphi - New diffeomorphic transform
%   ipsi - New complete affine+diffeomorphic mapping
%   bb   - New Bounding box
%   pf   - New pushed image
%   c    - New pushed voxel count
%   wmu  - New pulled template
%__________________________________________________________________________

    % --- Parse inputs
    p = inputParser;
    p.FunctionName = 'lsLatent';
    p.addRequired('model',  @(X) isstruct(X) && isfield(X, 'name'));
    p.addRequired('dz',     @checkarray);
    p.addRequired('z0',     @checkarray);
    p.addRequired('v0',     @checkarray);
    p.addRequired('llm0',   @isscalar);
    p.addRequired('mu',     @checkarray);
    p.addRequired('f',      @checkarray);
    p.addParameter('llz0',     nan,    @isscalar);
    p.addParameter('regz',     []);
    p.addParameter('A',        eye(4), @(X) isnumeric(X) && issame(size(X), [4 4]));
    p.addParameter('Mf',       eye(4), @(X) isnumeric(X) && issame(size(X), [4 4]));
    p.addParameter('Mmu',      eye(4), @(X) isnumeric(X) && issame(size(X), [4 4]));
    p.addParameter('nit',      6,      @isscalar);
    p.addParameter('itgr',     nan,    @isscalar);
    p.addParameter('prm',      [0.0001 0.001 0.2 0.05 0.2], @(X) length(X) == 5);
    p.addParameter('bnd',      0,      @(X) isscalar(X) && isnumeric(X));
    p.addParameter('match',    'pull', @ischar);
    p.addParameter('geod',     0,      @isscalar);
    p.addParameter('llv0',     nan,    @isscalar);
    p.addParameter('par',      false,  @isscalar);
    p.addParameter('loop',     '',     @(X) ischar(X) && any(strcmpi(X, {'slice', 'component', 'none', ''})));
    p.addParameter('pf',  nan, @(X) isnumeric(X) || isa(X, 'file_array'));
    p.addParameter('c',   nan, @(X) isnumeric(X) || isa(X, 'file_array'));
    p.addParameter('wa',  nan, @(X) isnumeric(X) || isa(X, 'file_array'));
    p.addParameter('wmu', nan, @(X) isnumeric(X) || isa(X, 'file_array'));
    p.addParameter('verbose',  false,  @isscalar);
    p.addParameter('debug',    false,  @isscalar);
    p.parse(model, dz, z0, v0, llm0, mu, f, varargin{:});
    llz0    = p.Results.llz0;
    llv0    = p.Results.llv0;
    geod    = p.Results.geod;
    regz    = p.Results.regz;
    A       = p.Results.A;
    Mf      = p.Results.Mf;
    Mmu     = p.Results.Mmu;
    nit     = p.Results.nit;
    itgr    = p.Results.itgr;
    prm     = p.Results.prm;
    bnd     = p.Results.bnd;
    par     = p.Results.par;
    loop    = p.Results.loop;
    verbose = p.Results.verbose;
    debug   = p.Results.debug;
    
    matchmode = p.Results.match;
    pf      = p.Results.pf;
    c       = p.Results.c;
    wa      = p.Results.wa;
    wmu     = p.Results.wmu;
    
    if debug, fprintf('* lsLatent\n'); end
    
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
    
    % --- Template voxel size
    vsmu = sqrt(sum(Mmu(1:3,1:3).^2)); 
    
    
    % --- Load some data (in case it is on disk)
    z0 = numeric(z0);
    v0 = numeric(v0);
    dz = numeric(dz);
    dv = reconstructVelocity('latent', dz, 'subspace', W, ...
        'debug', debug, 'par', par, 'loop', loop);
    
    % --- Set some default parameter value
    if isempty(regz)
        regz = precisionZ(W, vsmu, prm, bnd, 'debug', debug);
    end
    if isnan(llz0)
       llz0 = llPriorLatent(z0, regz, 'fast');
    end
    if isnan(llv0)
        if geod
            llv0 = geod * llPriorVelocity(v0, 'fast', 'vs', vsmu, 'prm', prm, 'bnd', bnd);
        else
            llv0 = 0;
        end
    end
    
    % --- Initialise line search
    armijo = 1;           % Armijo factor
    ll0    = llm0 + llz0 + llv0; % Log-likelihood (only parts that depends on z)
    dimf   = [size(f) 1 1];
    latf   = dimf(1:3);
    dimmu  = [size(mu) 1 1];
    latmu  = dimmu(1:3);
    
    if verbose
        printInfo('header');
        printInfo('initial', ll0, llm0, llz0, llv0);
    end
    
    % //!\\ Here, the pushed image is temporarily kept in memory. This
    % might use up a large ammount of RAM when a lot of classes are used.
    % I might need to store it in a temporary file_array.
    
    % --- Loop
    for i=1:nit
        z = single(z0 + dz / armijo);
        v = single(v0 + dv / armijo);
        iphi = exponentiateVelocity(v, 'iphi', 'itgr', itgr, 'vs', vsmu, 'prm', prm, 'bnd', bnd, 'debug', debug);
        ipsi = reconstructIPsi(A, iphi, 'lat', latf, 'Mf', Mf, 'Mmu', Mmu, 'debug', debug);
        if strcmpi(matchmode, 'push')
            [pf, c, bb] = pushImage(ipsi, f, latmu, 'par', par, 'loop', loop, 'debug', debug, 'output', {pf, c});
            llm = llMatching(model, mu, pf, c, 'bb', bb, 'par', par, 'loop', loop, 'debug', debug);
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
            llm = llMatching(model, wmu, f, 'par', par, 'loop', loop, 'debug', debug);
        end
        if geod
            llv = geod * llPriorVelocity(v, 'fast', 'vs', vsmu, 'prm', prm, 'bnd', bnd);
        else
            llv = 0;
        end
        llz = llPriorLatent(z, regz, 'fast', 'debug', debug);
        ll  = llm + llz + llv;
        
        if verbose, printInfo(i, ll0, llm, llz, llv); end
        
        if ll <= ll0
            if verbose, printInfo('failed'); end
            armijo = armijo * 2;
        else
            if verbose, printInfo('success'); end
            result     = struct;
            result.ok  = true;
            result.z   = z;
            result.v   = v;
            result.llz = llPriorLatent(z, regz, 'debug', debug);
            result.llm = llm;
            result.iphi = iphi;
            result.ipsi = ipsi;
            if strcmpi(matchmode, 'pull')
                [pf, c, bb] = pushImage(ipsi, f, latmu, 'par', par, 'loop', loop, 'debug', debug, 'output', {pf, c});
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

function printInfo(which, oll, llm, llz, llv)
    if ischar(which) 
        if strcmpi(which, 'header')
            fprintf('Z - LineSearch | Armijo  | %12s = %12s + %12s + %12s | %12s\n', 'RLL', 'LL-Match', 'RLL-Prior', 'RLL-Geod', 'LL-Diff');
        elseif strcmpi(which, 'initial')
            fprintf('Z - LineSearch | Initial | %12.6f = %12.6f + %12.6f + %12.6f \n', oll, llm, llz, llv);
        elseif strcmpi(which, 'failed')
            fprintf('| Failed\n');
        elseif strcmpi(which, 'success')
            fprintf('| Success\n');
        elseif strcmpi(which, 'end')
            fprintf('Z - LineSearch | Complete failure\n');
        end
    else
        fprintf('Z - LineSearch | Try %3d | %12.6f = %12.6f + %12.6f + %12.6f | %12.6f ', which, llm+llz+llv, llm, llz, llv, llm+llz+llv-oll);
    end
end
        