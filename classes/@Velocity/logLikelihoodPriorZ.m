function ll = logLikelihoodPriorZ(obj, z, regz, fast)
% FORMAT ll = obj.logLikelihoodPriorZ((z), (regz), (fast))
% (z)    - Latent coordinates [default: obj.z]
% (regz) - Inverse covariance matrix of the prior [default: obj.regz]
% (ll)   - Log-likelihood of the input coordinates
%          [if no argout: write to obj.llz]
%
% Returns log p(Z | W)
    
    % -- Default arguments
    if nargin < 4
        fast = false;
        if nargin < 3
            obj.computeRegZ();
            regz = obj.regz;
            if nargin < 2
                z = obj.z;
            end
        end
    end
    
    % --- Check up-to-date
    if nargout == 0 && obj.utd.llz && ~fast
        ll = obj.llz;
        return
    end
    
    % --- Check all arrays are ready to be used
    if ~obj.checkarray(z) || ~obj.checkarray(regz)
        ll = 0;
        if obj.Debug
            warning('Cannot compute log-likelihood: missing arrays\n');
        end
        return
    end
    if obj.Debug, fprintf('* logLikelihoodPriorZ\n'); end;
    
    % --- Load in memory
    z    = double(numeric(z));
    regz = double(numeric(regz));
    
    % --- Compute log-likelihood
    ll = z' * regz * z;
    
    % --- Add constants w.r.t. Z
    if ~fast
        K = size(z,1);    % Number of principal components
        ll = ll + K * log(2 * pi) ...
                - logDet(regz);
    end
    
    % --- Common factor
    ll = -0.5 * ll;
    
    if nargout == 0
        obj.llz = ll;
        obj.statusChanged('llz');
    end

end
