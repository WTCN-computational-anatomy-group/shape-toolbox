function ll = llPriorAffine(q, regq, varargin)
% FORMAT ll = llPriorAffine(q, regq, ('fast'))
%
% ** Required **
% q    - Affine parameters in the Lie-algebra.
% regq - Precision matrix of the Affine Lie-algebra parameters.
% ** Optional **
% fast - If 'fast', only compute terms that depend on q [default: false]
% ** Output **
% ll   - Log-likelihood of the affine parameters (prior term)
%
% Returns log p(q).
% Non-regularised affine parameters are supposed stored at the front
% of the q vector.
    

    % --- Parse inputs
    p = inputParser;
    p.FunctionName = 'llPriorAffine';
    p.addRequired('q',     @checkarray);
    p.addRequired('regq',  @checkarray);
    p.addOptional('fast',   '', @(X) ischar(X) && strcmpi(X, 'fast'));
    p.addParameter('debug',  false, @isscalar);
    p.parse(q, regq, varargin{:});
    fast  = p.Results.fast;
    debug = p.Results.debug;
    
    if debug, fprintf('* llPriorAffine\n'); end;
    
    
    % --- Load in memory
    q    = double(numeric(q));
    regq = double(numeric(regq));
    sub  = (size(q, 1) - size(regq, 1) + 1):size(q, 1);
    
    % --- Compute log-likelihood
    ll = q(sub)' * regq * q(sub);
    
    % --- Add constants w.r.t. Z
    if ~strcmpi(fast, 'fast')
        Q = size(regq,1);    % Number of affine parameters with prior
        ll = ll + Q * log(2 * pi) ...
                - spm_matcomp('LogDet', regq);
    end
    
    % --- Common factor
    ll = -0.5 * ll;

end
