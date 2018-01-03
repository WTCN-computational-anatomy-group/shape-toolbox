function ll = llPriorVelocity(v, varargin)
%__________________________________________________________________________
% Prior log-likelihood of the initial velocity 
%    log p(V | L)
%--------------------------------------------------------------------------
% FORMAT ll = llPriorVelocity(v, ('fast'), ...)
%
% REQUIRED
% --------
% v      - Initial velocity
%
% OPTIONAL
% --------
% fast   - If 'fast', only compute terms that depend on v [default: false]
%
% KEYWORD ARGUMENTS
% -----------------
% vs     - Voxel size of the velocity lattice [1 1 1]
% prm    - Parameters of the L operator (see spm_diffeo)
%          [0.0001 0.001 0.2 0.05 0.2]
% bnd    - L differential operator boundary conditions (0/1/2/3) [0]
% logdet - Log determinant of the kernel (so that it is not computed
%          several times)
%
% OUTPUT
% ------
% ll     - Log-likelihood of the initial velocity (prior term)
%__________________________________________________________________________
    

    % --- Parse inputs
    p = inputParser;
    p.FunctionName = 'llPriorVelocity';
    p.addRequired('v',       @checkarray);
    p.addOptional('fast',    '', @(X) ischar(X) && strcmpi(X, 'fast'));
    p.addParameter('vs',     [1 1 1], @(X) isscalar(X) || checkarray(X));
    p.addParameter('prm',    [0.0001 0.001 0.2 0.05 0.2], @checkarray);
    p.addParameter('bnd',       0, @(X) isscalar(X) && isnumeric(X));
    p.addParameter('logdet',  [], @isnumeric);
    p.addParameter('debug',  false, @isscalar);
    p.parse(v, varargin{:});
    fast = strcmpi(p.Results.fast, 'fast');
    vs    = p.Results.vs;
    prm   = p.Results.prm;
    bnd   = p.Results.bnd;
    ld    = p.Results.logdet;
    debug = p.Results.debug;
    
    if debug, fprintf('* llPriorVelocity\n'); end
    
    spm_diffeo('boundary', bnd);
    
    % --- Load in memory
    v = single(numeric(v));
    dim = [size(v) 1 1];
    lat = dim(1:3);
    
    % --- Discard bad voxels (they should not exist but just in case)
    if ~fast, count = sumall(isfinite(v)); end
    v(~isfinite(v)) = 0;
    
    % --- Compute log-likelihood
    m = spm_diffeo('vel2mom', v, double([vs prm]));
    ll = v(:)' * m(:);
    
    % --- Add constants w.r.t. R
    if ~fast
        ll = ll + count * log(2 * pi);
        
        if ~isempty(ld)
            [~, ld] = spm_shoot_greens('kernel', double(lat), double([vs prm]));
            ld = ld(1);
        end
        ll = ll - ld;
    end
    
    % --- Common factor
    ll = - 0.5 * ll;
    
end
