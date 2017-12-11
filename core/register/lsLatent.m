function [ok, z, llm, llz, v, iphi, pf, c, bb, ipsi] = lsLatent(model, dz, z0, v0, llm0, W, mu, f, varargin)
% FORMAT [ok, z, llm, llz, v, iphi, pf, c, bb, ipsi] = lsLatent(model, dz, z0, v0, llm0, W, mu, f, ...)
%
% ** Required **
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
% ** Keyword arguments **
% llz0 - Previous log-likelihood (prior term) [compute]
% regz - Precision matrix of the latent parameters [none]
% A    - Affine transform [eye(4)]
% Mf   - Image voxel-to-world mappinf [eye(4)]
% Mmu  - Template voxel-to-world mapping [eye(4)]
% nit  - Number of line-search iterations [6]
% itgr - Number of integration steps for geodesic shooting [auto]
% prm  - Differential operator parameters [0.0001 0.001 0.2 0.05 0.2]
% loop - How to split processing [auto]
% par  - If true, parallelise processing [false]
% ** Output **
% ok   - True if a better parameter value was found
% z    - New parameter value
% llm  - New log-likelihood (matching term)
% llz  - New log-likelihood (prior term)
% iphi - New diffeomorphic transform
% pf   - New pushed image
% c    - New pushed voxel count
% bb   - New Bounding box
% ipsi - New complete affine+diffeomorphic mapping
%
% Performs a line search along a direction to find better latent
% coordinates. The line search direction is usually found by Gauss-Newton.

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
    p.addParameter('llz0',     nan,    @iscalar);
    p.addParameter('regz',     []);
    p.addParameter('A',        eye(4), @(X) isnumeric(X) && issame(size(X), [4 4]));
    p.addParameter('Mf',       eye(4), @(X) isnumeric(X) && issame(size(X), [4 4]));
    p.addParameter('Mmu',      eye(4), @(X) isnumeric(X) && issame(size(X), [4 4]));
    p.addParameter('nit',      6,      @isscalar);
    p.addParameter('itgr',     nan,    @isscalar);
    p.addParameter('prm',      [0.0001 0.001 0.2 0.05 0.2], @(X) length(X) == 5);
    p.addParameter('par',      false,  @isscalar);
    p.addParameter('loop',     '',     @(X) ischar(X) && any(strcmpi(X, {'slice', 'component', 'none', ''})));
    p.addParameter('output',   []);
    p.addParameter('verbose',  false,  @isscalar);
    p.addParameter('debug',    false,  @isscalar);
    p.parse(model, dz, z0, v0, llm0, mu, f, varargin{:});
    llz0    = p.Results.llz0;
    regz    = p.Results.regz;
    A       = p.Results.A;
    Mf      = p.Results.Mf;
    Mmu     = p.Results.Mmu;
    nit     = p.Results.nit;
    itgr    = p.Results.itgr;
    prm     = p.Results.prm;
    par     = p.Results.par;
    loop    = p.Results.loop;
    output  = p.Results.output;
    verbose = p.Results.verbose;
    debug   = p.Results.debug;
    
    if debug, fprintf('* lsLatent\n'); end;
    
    % --- Template voxel size
    vsmu = sqrt(sum(Mmu(1:3,1:3).^2)); 
    
    % --- Set some default parameter value
    if isempty(regz)
        regz = precisionZ(W, vsmu, prm, 'debug', debug);
    end
    if isnan(llz0)
       llz0 = llPriorLatent(z0, regz, 'fast');
    end
    
    % --- Load some data (in case it is on disk)
    z0 = numeric(z0);
    dz = numeric(dz);
    dv = reconstructVelocity('latent', dz, 'subspace', W, ...
        'debug', debug, 'par', par, 'loop', loop);
    
    % --- Initialise line search
    armijo = 1;           % Armijo factor
    ok     = false;       % Found a better ll ?
    ll0    = llm0 + llz0; % Log-likelihood (only parts that depends on z)
    dimf   = [size(f) 1 1];
    latf   = dimf(1:3);
    dimmu  = [size(mu) 1 1];
    latmu  = dimmu(1:3);
    
    if verbose
        printInfo('header');
        printInfo('initial', ll0, llm0, llz0);
    end
    
    % //!\\ Here, the pushed image is temporarily kept in memory. This
    % might use up a large ammount of RAM when a lot of classes are used.
    % I might need to store it in a temporary file_array.
    
    % --- Loop
    for i=1:nit
        z = single(z0 + dz / armijo);
        v = single(numeric(v0) + dv / armijo);
        iphi = exponentiateVelocity(v, 'iphi', 'itgr', itgr, 'vs', vsmu, 'prm', prm, 'debug', debug);
        ipsi = reconstructIPsi(A, iphi, 'lat', latf, 'Mf', Mf, 'Mmu', Mmu, 'debug', debug);
        [pf, c, bb] = pushImage(ipsi, f, latmu, 'par', par, 'loop', loop, 'debug', debug);
        llm = llMatching(model, mu, pf, c, 'bb', bb, 'par', par, 'loop', loop, 'debug', debug);
        llz = llPriorLatent(z, regz, 'fast', 'debug', debug);
        ll  = llm + llz;
        
        if verbose, printInfo(i, ll0, llm, llz); end;
        
        if ll <= ll0
            if verbose, printInfo('failed'); end;
            armijo = armijo * 2;
        else
            if verbose, printInfo('success'); end;
            llz = llPriorLatent(z, regz, 'debug', debug);
            ok  = true;
            return
        end
    end
    
    if verbose, printInfo('end'); end;
    z   = z0;

end

function printInfo(which, oll, llm, llz)
    if ischar(which) 
        if strcmpi(which, 'header')
            fprintf('Z - LineSearch | Armijo  | %12s = %12s + %12s | %12s\n', 'RLL', 'LL-Match', 'RLL-Prior', 'LL-Diff');
        elseif strcmpi(which, 'initial')
            fprintf('Z - LineSearch | Initial | %12.6f = %12.6f + %12.6f \n', oll, llm, llz);
        elseif strcmpi(which, 'failed')
            fprintf('| Failed\n');
        elseif strcmpi(which, 'success')
            fprintf('| Success\n');
        elseif strcmpi(which, 'end')
            fprintf('Z - LineSearch | Complete failure\n');
        end
    else
        fprintf('Z - LineSearch | Try %3d | %12.6f = %12.6f + %12.6f | %12.6f ', which, llm+llz, llm, llz, llm+llz-oll);
    end
end
        