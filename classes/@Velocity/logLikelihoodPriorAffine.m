function ll = logLikelihoodPriorAffine(obj, q, regq, fast)
% FORMAT ll = obj.logLikelihoodPriorZ((z), (regz), (fast))
% (q)    - Latent coordinates [default: obj.q]
% (regq) - Inverse covariance matrix of the prior [default: obj.regq]
% (ll)   - Log-likelihood of the input coordinates
%          [if no argout: write to obj.llq]
%
% Returns log p(q)
    
    % -- Default arguments
    if nargin < 4
        fast = false;
        if nargin < 3
            regq = obj.RegAffine;
            if nargin < 2
                q = obj.q;
            end
        end
    end
    
    % --- Check up-to-date
    if nargout == 0 && obj.utd.llq && ~fast
        ll = obj.llq;
        return
    end
    
    % --- Check all arrays are ready to be used
    if ~obj.checkarray(q) || ~obj.checkarray(regq)
        ll = 0;
        if obj.Debug
            warning('Cannot compute log-likelihood: missing arrays\n');
        end
        return
    end
    if obj.Debug, fprintf('* logLikelihoodPriorAffine\n'); end;
    
    % --- Load in memory
    q    = double(numeric(q));
    regq = double(numeric(regq));
    sub  = (size(q, 1) - size(regq, 1) + 1):size(q, 1);
    
    % --- Compute log-likelihood
    ll = q(sub)' * regq * q(sub);
    
    % --- Add constants w.r.t. Z
    if ~fast
        Q = size(regq,1);    % Number of affine parameters with prior
        ll = ll + Q * log(2 * pi) ...
                - logDet(regq);
    end
    
    % --- Common factor
    ll = -0.5 * ll;
    
    if nargout == 0
        obj.llq = ll;
        obj.statusChanged('llq');
    end

end
