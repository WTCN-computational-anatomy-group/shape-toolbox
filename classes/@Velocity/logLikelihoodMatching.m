function ll = logLikelihoodMatching(obj, varargin)
% FORMAT (ll) = obj.logLikelihoodMatching('res', (res))
% res    - Residuals [default: use obj.res]
% (ll)   - Log-likelihood [if no argout: write to obj.llm]
% > Compute from already computed residuals
%
% FORMAT (ll) = obj.logLikelihoodMatching('res', varargin)
% args   - See following formats
% (ll)   - Log-likelihood [if no argout: write to obj.llm]
% > Compute complete residuals based on the given arguments, then compute
% log-likelihood from residuals. If 'res' is not used, more memory/space
% efficient algorithms -- which do not require the complete residuals to
% exist -- will be used.
%
% FORMAT (ll) = obj.logLikelihoodMatching(mu, f)
% f    - Observed image aligned with mu
% mu   - Template image aligned with f
% (ll) - Log-likelihood [if no argout: write to obj.llm]
% > Warping is done before call, loglikelihood is computed as if in 
% native space.
%
% FORMAT (ll) = obj.logLikelihoodMatching(ipsi, mu, f, (type))
% ipsi   - Inverse transform (warps mu to f)
% f      - Native observed image
% mu     - Native template
% type   - 'warp' > Template is warped in image space
%          'push' > Image is pushed to template space
% (ll)   - Log-likelihood [if no argout: write to obj.llm]
%
% FORMAT (ll) = obj.logLikelihoodMatching(ipsi, (type))
% ipsi   - Inverse transform (warps mu to f)
% type   - 'warp' > Template is warped in image space
%          'push' > Image is pushed to template space
% (ll)   - Log-likelihood [if no argout: write to obj.llm]
% > Use obj.mu and obj.f
%
% FORMAT (ll) = obj.logLikelihoodMatching((type))
% type   - 'warp' > Template is warped in image space
%          'push' > Image is pushed to template space
% (ll)   - Log-likelihood [if no argout: write to obj.llm]
% > Use (obj.wmu/obj.f or obj.mu/obj.pf/obj.pvox)
%
% obj   - Velocity object.
%         * Use property 'MatchingTerm'.
%         * May use properties 'Interpolation', 'Normal.s', 'Laplace.b'.

    if nargout == 0 && obj.utd.llm
        ll = obj.llm;
        return
    end
    
    use_res = false;
    if nargin == 0
        use_res = true;
    end
    if numel(varargin) > 0 && ischar(varargin{1}) ...
            && strcmpi(varargin{1}, 'res')
        varargin = varargin(2:end);
        use_res = true;
    end
    
    % --------------------
    % - If use residuals -
    % --------------------
    if use_res
        
        % --- Provided residuals
        if numel(varargin) == 1
            res = varargin{1};
        else
            % --- Read optional type argument
            if numel(varargin) > 0 && ischar(varargin{end})
                type = varargin{end};
            else
                type = 'warp';
            end
            % --- Compute residuals
            if nargout == 0
                obj.residuals(varargin{:});
                if strcmpi(type, 'warp')
                    res = obj.res;
                else
                    res = obj.pres;
                end
            else
                res = obj.residuals(varargin{:});
            end
        end
        ll = logLikelihoodResiduals(res);
    
    % ----------------------------------
    % - Else (used during line search) -
    % ----------------------------------
    else
        % --- Read optional type argument
        if numel(varargin) > 0 && ischar(varargin{end})
            type = varargin{end};
            varargin = varargin(1:end-1);
        else
            type = 'warp';
        end
    
        % --- Select case
        switch numel(varargin)
            case 0
            % ll = obj.logLikelihood((type))
                if strcmp(type, 'warp')
                    obj.reconstructWarpedTemplate();
                    ll = logLikelihoodWarped(obj, obj.wmu, obj.f);
                else
                    obj.reconstructTemplate();
                    obj.pushImage();
                    ll = logLikelihoodPushed(obj, obj.mu, obj.pf, obj.pvox);
                end

            case 1
            % ll = obj.logLikelihood(iphi, (type))
            % ipsi = varargin{1}
                if strcmp(type, 'warp')
                    wmu = obj.warpTemplate(varargin{1});
                    wmu = obj.reconstructWarpedTemplate(wmu);
                    ll = logLikelihoodWarped(obj, wmu, obj.f);
                else
                    obj.reconstructTemplate();
                    [pf, pvox] = obj.pushImage(varargin{1});
                    ll = logLikelihoodPushed(obj, obj.mu, pf, pvox);
                end

            case 2
            % ll = obj.logLikelihood(mu, f)
            % mu = varargin{1}
            % f  = varargin{2}
                ll = logLikelihoodWarped(obj, varargin{1}, varargin{2});

            case 3
            % ll = obj.logLikelihood(iphi, mu, f, (type))
            % ipsi = varargin{1}
            % mu   = varargin{2}
            % f    = varargin{3}
                if strcmp(type, 'warp')
                    wmu = obj.warpTemplate(varargin{1}, varargin{2});
                    wmu = obj.reconstructWarpedTemplate(wmu);
                    ll = logLikelihoodWarped(obj, wmu, varargin{3});
                else
                    [pf, pvox] = obj.pushImage(varargin{3}, varargin{1});
                    ll = logLikelihoodPushed(obj, varargin{2}, pf, pvox);
                end

            otherwise
                error('Too many arguments');
        end
    end
    
    % --- Write to disk
    if nargout == 0
        obj.llm     = ll;
        obj.utd.llm = true;
        obj.statusChanged('llm');
    end
end
% =========================================================================
function ll = logLikelihoodWarped(obj, mu, f)
% FORMAT ll = logLikelihoodWarped(obj, mu, f)
% obj   - Velocity object
%         Used properties: MatchingTerm, (Normal.s), (Laplace.b)
% mu    - Warped template in image space
% f     - Native image
% ll    - Log-likelihood

    % --- Check if input arrays are usable
    if ~obj.checkarray(mu) || ~obj.checkarray(f)
        ll = 0;
        if obj.Debug
            warning('Cannot compute matching log-likelihood (warped): missing arrays');
        end
        return
    end
    if obj.Debug, fprintf('* logLikelihoodWarped\n'); end;
    
    % --- Read dim info
    dim = [size(mu) 1 1 1 1];
    dim = dim(1:4);
    nc  = dim(4);
    
    % --- Deal with different matching terms
    switch lower(obj.MatchingTerm)
        case {'normal', 'gaussian', 'l2'}
            % E = log p(f | mu)
            %   = - sum{k,i} [1/(2s(k))  (mu(k)(xi) - f(k)(xi))^2 + Cte]
            
            s = obj.Normal.s;   % < sigma^2
            if length(s) == 1
                s = s * ones(1, nc);
            end
            ll = 0;
            for k=1:nc
                % Load one class (if on disk, saves memory)
                mu1 = single(mu(:,:,:,k));
                f1  = single(f(:,:,:,k));
                msk = isfinite(f1) & isfinite(mu1);
                count = sumall(msk > 0);
                % Compute log-likelihood
                ll = ll ...
                    - 0.5 * count * log(s(k)) ...
                    - 0.5 * count * log(2 * pi) ...
                    - 0.5 / s(k) * sumall((mu1(msk) - f1(msk)).^2);
                
                clear f1 mu1 msk
            end
            clear f mu
            
        case {'laplace', 'l1'}
            % E = log p(f | mu)
            %   = - sum{k,i} [1/b(k) |mu(k)(xi) - f(k)(xi)| + Cte]
    
            b = obj.Laplace.b;
            if length(b) == 1
                b = b * ones(1, nc);
            end
            ll = 0;
            for k=1:nc
                % Load one class (if on disk, saves memory)
                mu1 = single(mu(:,:,:,k));
                f1  = single(f(:,:,:,k));
                msk = isfinite(f1) & isfinite(mu1);
                count = sumall(msk > 0);
                % Compute log-likelihood
                ll = ll ...
                    - count * log(2 * b(k)) ...
                    - 1/b * sumall(abs(mu1(msk) - f1(msk)));
                clear f1 mu1 msk
            end
            clear f mu
            
        case {'bernoulli', 'binomial', 'binary'}
            % Only one class !
            
            % Here, the input template is actually encoded in some kind of  
            % log space (that we will note a). Actual probabilities mu are 
            % mapped back with a sigmoid function:
            % mu = exp(a) / (1 + exp(a))
            
            % E = log p(f | mu, v)
            %   = sum{i} [f(xi) log(mu(xi)) + (1 - f(xi)) log(1 - mu(xi))]
            
            mu = single(mu(:,:,:,1));
            f  = single(f(:,:,:,1));
            msk = isfinite(f) & isfinite(mu);
            % Just in case: correct for interpolation
            mu = max(eps('single'), min(1-eps('single'), mu));
            % Compute log-likelihood
            ll = sumall(f(msk) .* log(mu(msk)) + (1 - f(msk)) .* log(1 - mu(msk)));
            clear f mu msk
            
        case {'categorical', 'multinomial'}
            % Here, the input template is actually encoded in some kind of  
            % log space (that we will note a). Actual probabilities mu are 
            % mapped back with a softmax function:
            % mu[k] = exp(a[k]) / (sum{l} exp(a[l]))

            % E = log p(f | mu, v)
            %   = sum{i,k} [ f(k)(xi) log(mu(k)(xi)) ]
            
            for k=1:nc
                mu1 = single(mu(:,:,:,k));
                f1  = single(f(:,:,:,k));
                msk = isfinite(f1) & isfinite(mu1);
                mu1 = max(eps('single'), mu1);  % Correct for interpolation
                ll = ll + sumall(f1(msk) .* log(mu1(msk)));
                clear f1 mu1 msk
            end
            clear f mu
    end
end
% =========================================================================
function ll = logLikelihoodPushed(obj, mu, f, c)
% FORMAT ll = logLikelihoodWarped(obj, mu, f, c)
% obj   - Velocity object
%         Used properties: MatchingTerm, (Normal.s), (Laplace.b)
% mu    - Template image in templatespace
% f     - Observed image pushed in template space
% c     - Pushed voxel count
% ll    - Log-likelihood

    % --- Check if input arrays are usable
    if ~obj.checkarray(mu) || ~obj.checkarray(f) || ~obj.checkarray(c)
        ll = 0;
        if obj.Debug
            warning('Cannot compute matching log-likelihood (pushed): missing arrays\n');
        end
        return
    end
    if obj.Debug, fprintf('* loglikelihoodPushed\n'); end;
    
    % --- Read dim info
    dim = [size(mu) 1 1 1 1];
    dim = dim(1:4);
    % lat = dim(1:3);
    nc  = dim(4);
    
    % --- Deal with different matching terms
    switch lower(obj.MatchingTerm)
        case {'normal', 'gaussian', 'l2'}
            % E = log p(f | mu)
            %   = - sum{k,i} [1/(2s(k)) c(xi) (mu(k)(xi) - f(k)(xi)/c(xi))^2 + Cte]
    
            s = obj.Normal.s;   % < sigma^2
            if length(s) == 1
                s = s * ones(1, nc);
            end
            c = single(numeric(c)); % Load count (if on disk, saves i/o)
            ll = 0;
            for k=1:nc
                % Load one class (if on disk, saves memory)
                mu1 = single(mu(:,:,:,k));
                f1  = single(f(:,:,:,k));
                msk = isfinite(mu1) & isfinite(f1) & (c > 0);
                count = sumall(msk > 0);
                % Compute log-likelihood
                ll = ll ...
                    - 0.5 * count * log(s(k)) ...
                    - 0.5 * count * log(2 * pi) ...
                    - 0.5 /s(k) * sumall((mu1(msk) .* c(msk) - f1(msk)).^2);
                clear f1 mu1 msk
            end
            clear f mu c
            
        case {'laplace', 'l1'}
            % E = log p(f | mu)
            %   = - sum{k,i} [1/b(k) c(xi) |mu(k)(xi) - f(k)(xi)/c(xi)| + Cte]
    
            b = obj.Laplace.b;
            if length(b) == 1
                b = b * ones(1, nc);
            end
            c = single(numeric(c)); % Load count (if on disk, saves i/o)
            ll = 0;
            for k=1:nc
                % Load one class (if on disk, saves memory)
                mu1 = single(mu(:,:,:,k));
                f1  = single(f(:,:,:,k));
                msk = isfinite(mu1) & isfinite(f1) & (c > 0);
                count = sumall(msk > 0);
                % Compute log-likelihood
                ll = ll ...
                    - count * log(2 * b(k)) ...
                    - 1/b * sumall(abs(mu1(msk) .* c(msk) - f1(msk)));
                clear f1 mu1 msk
            end
            clear f mu c
            
        case {'bernoulli', 'binomial', 'binary'}
            % Only one class !
            
            % Here, the input template is actually encoded in some kind of  
            % log space (that we will note a). Actual probabilities mu are 
            % mapped back with a sigmoid function:
            % mu = exp(a) / (1 + exp(a))
            
            % E = log p(f | mu, v)
            %   = sum{i} c(xi) [(pf(xi)/c(xi)) log(mu(xi)) + (1 - pf(xi)/c(xi)) log(1 - mu(xi))]
            %   = sum{i} [pf(xi) log(mu(xi)) + (c(xi) - pf(xi)) log(1 - mu(xi))]
            
            c  = single(numeric(c));
            mu = single(mu(:,:,:,1));
            f  = single(f(:,:,:,1));
            msk = isfinite(mu) & isfinite(f) & (c > 0);
            mu = max(eps('single'), min(1-eps('single'), mu));
            % Compute log-likelihood
            ll = sumall(f(msk) .* log(mu(msk)) + (c(msk) - f(msk)) .* log(1 - mu(msk)));
            clear f mu c msk
            
        case {'categorical', 'multinomial'}
            % Here, the input template is actually encoded in some kind of  
            % log space (that we will note a). Actual probabilities mu are 
            % mapped back with a softmax function:
            % mu[k] = exp(a[k]) / (sum{l} exp(a[l]))
            
            % E = log p(f | mu, v)
            %   = sum{i,k} [ c(xi) (pf(k)(xi)/c(xi)) log(mu(k)(xi)) ]
            %   = sum{i,k} [ pf(k)(xi) log(mu(k)(xi)) ]
            ll = 0;
            for k=1:nc
                % Load one class (if on disk, saves memory)
                mu1 = single(mu(:,:,:,k));
                f1  = single(f(:,:,:,k));
                msk = isfinite(mu1) & isfinite(f1) & (c > 0);
                mu1 = max(eps('single'), mu1);
                ll = ll + sumall(f1(msk) .* log(mu1(msk)));
                clear f1 mu1 msk
            end
            clear f mu c
        otherwise
            error('Unknown distribution ''%s''', obj.MatchingTerm);
    end
end
% =========================================================================
function ll = logLikelihoodResiduals(obj, res)
% res - Residuals
% ll  - Log-likelihood (matching term)
% obj - Velocity object. This function will use properties 'MatchingTerm',
% 'Normal.s', 'Laplace.b'
%
% Compute log-likelihood from the image of residuals

    % 'normal':      res(k) = 0.5/s(k) .* (f(k) - mu(k)).^2
    % 'laplace':     res(k) = 1/b(k) .*abs(f(k) - mu(k))
    % 'binomial':    res = f .* mu + (1 - f) .* log(1 - exp(mu));
    % 'multinomial': res = sum(k)[f(k) .* mu(k)] - log(sum(k)[exp(mu(k))])
    
    if ~obj.checkarray(res)
        ll = 0;
        if obj.Debug
            warning('Cannot compute matching log-likelihood from residuals. Missing arrays')
        end
        return
    end
    if obj.Debug, fprintf('* logLikelihoodResiduals\n'); end;
    
    % --- Dim info
    nc = size(res, 4);
    
    % --- Select case
    switch lower(obj.MatchingTerm)
        case {'normal', 'gaussian', 'l2'}
            s = obj.Normal.s;   % < sigma^2
            if length(s) == 1
                s = s * ones(1, nc);
            end
            ll = 0;
            for k=1:nc
                % Load one class (if on disk, saves memory)
                res1  = single(res(:,:,:,k));
                mask  = isfinite(res1);
                count = sumall(mask > 0);
                % Compute log-likelihood
                ll = ll ...
                    - 0.5 * count * log(s(k)) ...
                    - 0.5 * count * log(2 * pi) ...
                    - 0.5 * sumall(res1(mask));
                clear res1
            end
            
        case {'laplace', 'l1'}
            b = obj.Laplace.b;
            if length(b) == 1
                b = b * ones(1, nc);
            end
            ll = 0;
            for k=1:nc
                % Load one class (if on disk, saves memory)
                res1  = single(res(:,:,:,k));
                mask  = isfinite(res1);
                count = sumall(mask > 0);
                % Compute log-likelihood
                ll = ll ...
                    - count * log(2 * b(k)) ...
                    - sumall(res1(mask));
                clear res1
            end
            
        case {'bernoulli', 'binomial', 'binary', 'categorical', 'multinomial'}
            ll = 0;
            for k=1:nc
                % Load one class (if on disk, saves memory)
                res1 = single(res(:,:,:,k));
                mask = isfinite(res1);
                % Compute log-likelihood
                ll = ll + sumall(res1(mask));
                clear res1
            end
    end
end