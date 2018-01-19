function ll = llPriorLatent(z, regz, varargin)
% FORMAT ll = llPriorLatent(z, regz, ('fast'))
%
% ** Required **
% z    - Latent coordinates
% regz - Precision matrix of the latent space (W'*L*W)
% ** Optional **
% fast - If 'fast', only compute terms that depend on z [default: false]
% ** Output **
% ll   - Log-likelihood of the latent coordinates (prior term)
%
% Returns log p(Z | W)
    

    % --- Parse inputs
    p = inputParser;
    p.FunctionName = 'llPriorLatent';
    p.addRequired('z',     @checkarray);
    p.addRequired('regz',  @checkarray);
    p.addOptional('fast',   '', @(X) ischar(X) && strcmpi(X, 'fast'));
    p.addParameter('debug',  false, @isscalar);
    p.parse(z, regz, varargin{:});
    
    if p.Results.debug, fprintf('* llPriorLatent\n'); end;
    
    % --- Load in memory
    z    = double(numeric(z));
    regz = double(numeric(regz));
    
    % --- Compute log-likelihood
    ll = z' * regz * z;
    
    % --- Add constants w.r.t. Z
    if ~strcmpi(p.Results.fast, 'fast')
        K = size(z,1);    % Number of principal components
        ll = ll + K * log(2 * pi) - spm_matcomp('LogDet', regz);
    end
    
    % --- Common factor
    ll = -0.5 * ll;

end
