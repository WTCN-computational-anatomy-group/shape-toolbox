function ll = llPriorVelocity(v, varargin)
% FORMAT ll = llPriorVelocity(v, ('fast'), ('vs',   vs), ('prm',  prm))
%
% ** Required **
% v    - Initial velocity
% ** Optional **
% fast - If 'fast', only compute terms that depend on v [default: false]
% ** Keyword arguments **
% vs   - Voxel size of the velocity lattice [1 1 1]
% prm  - Parameters of the L operator (see spm_diffeo)
%        [0.0001 0.001 0.2 0.05 0.2]
% ** Output **
% ll   - Log-likelihood of the initial velocity (prior term)
%
% Returns log p(V | L)
    

    % --- Parse inputs
    p = inputParser;
    p.FunctionName = 'llPriorVelocity';
    p.addRequired('v',       @checkarray);
    p.addOptional('fast',    '', @(X) ischar(X) && strcmpi(X, 'fast'));
    p.addParameter('vs',     [1 1 1], @(X) isscalar(X) || checkarray(X));
    p.addParameter('prm',    [0.0001 0.001 0.2 0.05 0.2], @checkarray);
    p.addParameter('debug',  false, @isscalar);
    p.parse(v, varargin{:});
    fast = strcmpi(p.Results.fast, 'fast');
    
    if p.Results.debug, fprintf('* llPriorVelocity\n'); end;
    
    % --- Load in memory
    v = single(numeric(v));
    dim = [size(v) 1 1];
    lat = dim(1:3);
    
    % --- Discard bad voxels (they should not exist but just in case)
    if ~fast, count = sumall(isfinite(v)); end;
    v(~isfinite(v)) = 0;
    
    % --- Compute log-likelihood
    m = spm_diffeo('vel2mom', v, [p.Results.vs p.Results.prm]);
    ll = v(:)' * m(:);
    
    % --- Add constants w.r.t. R
    if ~fast
        ll = ll + count * log(2 * pi);
        
        [~, ld] = spm_shoot_greens('kernel', lat, [p.Results.vs p.Results.prm]);
        ll = ll - 0.5*ld(1);
%         if obj.Debug
%             warning(['LL: discard 0.5 * det(L). ' ...
%                      'Don''t know how to compute it (with the kernel?).']);
%         end;
    end
    
    % --- Common factor
    ll = - 0.5 * ll;
    
end
