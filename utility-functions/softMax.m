function [mu, ll] = softMax(a, f, fast, c)
% FORMAT [sig, ll] = softMax(a, (f), (fast), (c))
% a      - Log of (reconstructed) TPMs
% If log-likelihood needed:
% (f)    - Pushed image in the template space
% (fast) - *     1: Assume probabilities are normalized
%                   (sum_k f_k and sum_k exp(a_k) should be 1)
%          * [1 0]: Insure normalization of the TPMs (softmax)
%          * [0 1]: Insure normalization of the observed
%          *     0: Insure both normalizations (default)
% (c)    - Pushed voxel count
%          (only useful if fast is 1 and f was pushed to template space)
% Output:
% mu     - SoftMax of the TPMs (exp(a_k)/sum(exp(a_j)))
% (ll)   - Log-likelihood of the observed: log p(f | mu) = f log mu

    % --- Default options
    if nargin < 4
        c = 1;
    end
    if nargin < 3
        fast = [0 0];
    elseif numel(fast) == 1
        fast = padarray(fast, [0 1], 'replicate', 'post');
    end
    
    % --- Load from disk
    if isa(a, 'file_array')
        a = a(:,:,:,:);
    end
    if nargin > 1 && isa(f, 'file_array')
        f = f(:,:,:,:);
    end
    
    % --- ??
    mxa = max(a, [], 4);
    a   = bsxfun(@minus, a, mxa);
    
    % --- Softmax exponential
    mu = exp(a);
    if ~fast(1)
        s = sum(mu, 4);
    end
    
    % --- log-likelihood
    ll  = 0;
    if nargin >=2 && nargout >= 2
        f(~isfinite(f)) = 0;
        if fast(1)
            ll = sumall( f .* a );
        elseif fast(2)
            ll = sumall( sum(f .* a, 4) - c .* log(s) );
        else
            ll = sumall( sum(f .* a, 4) - log(s) .* sum(f, 4) );
        end
    end
    clear a f
    
    % --- Softmax normalization
    if ~fast(1)
        mu = bsxfun(@rdivides, mu, s);
    end
end