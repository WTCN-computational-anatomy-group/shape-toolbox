function [ok, r, llm, llr, iphi, pf, c, bb, ipsi, v] = lsVelocity(model, dr, r0, llm0, mu, f, varargin)
%__________________________________________________________________________
% 
% Performs a line search along a direction to find better residual initial
% velocity. The line search direction is usually found by Gauss-Newton.
%
%--------------------------------------------------------------------------
%
% FORMAT [ok, ...] = lsVelocity(model, dr, r0, llm0, mu, f, ...)
%
% REQUIRED
% --------
% model - Structure with fields:
%           * 'name'    : 'normal', 'laplace', 'bernoulli' or 'categorical'
%           * ('sigma2'): Normal variance  [1]
%           * ('b')     : Laplace variance [1]
% dr   - Line search direction (ascent if maximising, descent if minimising)
% r0   - Previous residual field
% llm0 - Previous log-likelihood (matching term)
% W    - Principal subspace
% mu   - Template (in native space)
% f    - Image (in native space)
%
% KEYWORD ARGUMENTS
% -----------------
% v0   - Previous initial velocity [r0]
% llz0 - Previous log-likelihood (prior term) [compute]
% sigma- Variance captured by the residual field [1]
% A    - Affine transform [eye(4)]
% Mf   - Image voxel-to-world mappinf [eye(4)]
% Mmu  - Template voxel-to-world mapping [eye(4)]
% nit  - Number of line-search iterations [6]
% itgr - Number of integration steps for geodesic shooting [auto]
% prm  - Differential operator parameters [0.0001 0.001 0.2 0.05 0.2]
% bnd  - Differential operator boundary conditions (0/1/2/3) [0]
% loop - How to split processing [auto]
% par  - If true, parallelise processing [false]
%
% OUTPUT
% ------
% ok   - True if a better parameter value was found
% v    - New initial velocity
% llm  - New log-likelihood (matching term)
% llr  - New log-likelihood (prior term)
% iphi - New diffeomorphic transform
% pf   - New pushed image
% c    - New pushed voxel count
% ipsi - New complete affine+diffeomorphic mapping
%
%--------------------------------------------------------------------------
% 
% To perform the usual diffeomorphic registration (without principal 
% geodesic), set sigma == 1 and v0 = r0. Then, we always have r = v.
%__________________________________________________________________________

    % --- Parse inputs
    p = inputParser;
    p.FunctionName = 'lsVelocity';
    p.addRequired('model',  @(X) isstruct(X) && isfield(X, 'name'));
    p.addRequired('dr',     @checkarray);
    p.addRequired('r0',     @checkarray);
    p.addRequired('llm0',   @isscalar);
    p.addRequired('mu',     @checkarray);
    p.addRequired('f',      @checkarray);
    p.addParameter('v0',       [],     @checkarray);
    p.addParameter('llr0',     nan,    @isscalar);
    p.addParameter('sigma',    1,      @isscalar);
    p.addParameter('A',        eye(4), @(X) isnumeric(X) && issame(size(X), [4 4]));
    p.addParameter('Mf',       eye(4), @(X) isnumeric(X) && issame(size(X), [4 4]));
    p.addParameter('Mmu',      eye(4), @(X) isnumeric(X) && issame(size(X), [4 4]));
    p.addParameter('nit',      6,      @isscalar);
    p.addParameter('itgr',     nan,    @isscalar);
    p.addParameter('prm',      [0.0001 0.001 0.2 0.05 0.2], @(X) length(X) == 5);
    p.addParameter('shoot',    [0.0001 0.001 0.2 0.05 0.2], @(X) length(X) == 5);
    p.addParameter('bnd',      0, @(X) isscalar(X) && isnumeric(X));
    p.addParameter('par',      false,  @isscalar);
    p.addParameter('loop',     '',     @(X) ischar(X) && any(strcmpi(X, {'slice', 'component', 'none', ''})));
    p.addParameter('output',   []);
    p.addParameter('verbose',  false,  @isscalar);
    p.addParameter('debug',    false,  @isscalar);
    p.parse(model, dr, r0, llm0, mu, f, varargin{:});
    v0      = p.Results.v0;
    llr0    = p.Results.llr0;
    sigma   = p.Results.sigma;
    A       = p.Results.A;
    Mf      = p.Results.Mf;
    Mmu     = p.Results.Mmu;
    nit     = p.Results.nit;
    itgr    = p.Results.itgr;
    prm     = p.Results.prm;
    shoot   = p.Results.shoot;
    bnd     = p.Results.bnd;
    par     = p.Results.par;
    loop    = p.Results.loop;
    verbose = p.Results.verbose;
    debug   = p.Results.debug;
    
    if debug, fprintf('* lsVelocity\n'); end
    
    % --- Template voxel size
    vsmu = sqrt(sum(Mmu(1:3,1:3).^2)); 
    
    % --- Set some default parameter value
    if isempty(v0)
        v0 = r0;
    end
    if isnan(llr0)
       llr0 = llPriorVelocity(r0, 'fast', 'vs', vsmu, 'prm', prm, 'bnd', bnd, 'debug', debug);
    end
    
    % --- Load some data (in case it is on disk)
    r0 = numeric(r0);
    v0 = numeric(v0);
    dr = numeric(dr);
    
    % --- Initialise line search
    armijo = 1;           % Armijo factor
    ok     = false;       % Found a better ll ?
    ll0    = llm0 + llr0; % Log-likelihood (only parts that depends on v)
    dimf   = [size(f) 1 1];
    latf   = dimf(1:3);
    dimmu  = [size(mu) 1 1];
    latmu  = dimmu(1:3);
    
    if verbose
        printInfo('header');
        printInfo('initial', ll0, llm0, llr0);
    end
    
    % //!\\ Here, the pushed image is temporarily kept in memory. This
    % might use up a large ammount of RAM when a lot of classes are used.
    % I might need to store it in a temporary file_array.
    
    % --- Loop
    for i=1:nit
        r = single(r0 + dr / armijo);
        v = single(v0 + sigma * dr / armijo);
        iphi = exponentiateVelocity(v, 'iphi', 'itgr', itgr, 'vs', vsmu, 'prm', shoot, 'bnd', bnd, 'debug', debug);
        ipsi = reconstructIPsi(A, iphi, 'lat', latf, 'Mf', Mf, 'Mmu', Mmu, 'debug', debug);
        [pf, c, bb] = pushImage(ipsi, f, latmu, 'par', par, 'loop', loop, 'debug', debug);
        llm = llMatching(model, mu, pf, c, 'bb', bb, 'par', par, 'loop', loop, 'debug', debug);
        llr = llPriorVelocity(r,  'fast', 'vs', vsmu, 'prm', prm, 'bnd', bnd, 'debug', debug);
        ll  = llm + llr;
        
        if verbose, printInfo(i, ll0, llm, llr); end
        
        if ll <= ll0
            if verbose, printInfo('failed'); end
            armijo = armijo * 2;
        else
            if verbose, printInfo('success'); end
            llr = llPriorVelocity(r, 'vs', vsmu, 'prm', prm, 'bnd', bnd, 'debug', debug);
            ok  = true;
            return
        end
    end
    
    if verbose, printInfo('end'); end
    r   = r0;
    v   = v0;

end

function printInfo(which, oll, llm, llr)
    if ischar(which) 
        if strcmpi(which, 'header')
            fprintf('R - LineSearch | Armijo  | %12s = %12s + %12s | %12s\n', 'RLL', 'LL-Match', 'RLL-Prior', 'LL-Diff');
        elseif strcmpi(which, 'initial')
            fprintf('R - LineSearch | Initial | %12.6f = %12.6f + %12.6f \n', oll, llm, llr);
        elseif strcmpi(which, 'failed')
            fprintf('| Failed\n');
        elseif strcmpi(which, 'success')
            fprintf('| Success\n');
        elseif strcmpi(which, 'end')
            fprintf('R - LineSearch | Complete failure\n');
        end
    else
        fprintf('R - LineSearch | Try %3d | %12.6f = %12.6f + %12.6f | %12.6f ', which, llm+llr, llm, llr, llm+llr-oll);
    end
end
        