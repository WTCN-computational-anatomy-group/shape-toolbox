function varargout = residuals(obj, varargin)
% FORMAT ([r, dr]) = obj.residuals(mu, f)
% f    - Observed image aligned with mu
% mu   - Template image aligned with f
% (r)  - Residuals used for matching term [if no argout: write to obj.res]
% (dr) - Residuals used for gradient of matching term
%        [if no argout: write to obj.dres]
% > Warping is done before call, loglikelihood is computed as if in 
% native space.
%
% FORMAT ([r, dr]) = obj.residuals(ipsi, mu, f, (type))
% ipsi   - Inverse transform (warps mu to f)
% f      - Native observed image
% mu     - Native template
% type   - 'warp' > Template is warped in image space
%          'push' > Image is pushed to template space
% (r)    - Residuals used for matching term
%          [if no argout: write to obj.res (warp) / obj.pres (push)]
% (dr)   - Residuals used for gradient of matching term
%          [if no argout: write to obj.presd (push)]
%
% FORMAT ([r, dr]) = obj.residuals(ipsi, (type))
% ipsi   - Inverse transform (warps mu to f)
% type   - 'warp' > Template is warped in image space
%          'push' > Image is pushed to template space
% (r)    - Residuals used for matching term
%          [if no argout: write to obj.res (warp) / obj.pres (push)]
% (dr)   - Residuals used for gradient of matching term
%          [if no argout: write to obj.presd (push)]
% > Use obj.mu and obj.f
%
% FORMAT ([r, dr]) = obj.residuals((type))
% type   - 'warp' > Template is warped in image space
%          'push' > Image is pushed to template space
% (r)    - Residuals used for matching term
%          [if no argout: write to obj.res (warp) / obj.pres (push)]
% (dr)   - Residuals used for gradient of matching term
%          [if no argout: write to obj.presd (push)]
% > Use (obj.wmu/obj.f or obj.mu/obj.pf/obj.pvox)
%
% obj   - Velocity object.
%         * Use property 'MatchingTerm'.
%         * May use properties 'Interpolation', 'Normal.s', 'Laplace.b'.

    % --- Read optional type argument
    if numel(varargin) > 0 && ischar(varargin{end})
        type = varargin{end};
        varargin = varargin(1:end-2);
    else
        type = 'warp';
    end
    
    % --- Select case
    switch numel(varargin)
        case 0
        % ll = obj.residuals((type))
            if strcmp(type, 'warp')
                obj.reconstructWarpedTemplate();
                [varargout{1:nargout}] = residualsWarped(obj, obj.wmu, obj.f);
            else
                obj.reconstructTemplate();
                obj.pushImage();
                [varargout{1:nargout}] = residualsPushed(obj, obj.mu, obj.pf, obj.pvox);
            end
            
        case 1
        % ll = obj.residuals(iphi, (type))
        % ipsi = varargin{1}
            if strcmp(type, 'warp')
                wmu = obj.warpTemplate(varargin{1});
                wmu = obj.reconstructWarpedTemplate(wmu);
                [varargout{1:nargout}] = residualsWarped(obj, wmu, obj.f);
            else
                obj.reconstructTemplate();
                [pf, pvox] = obj.pushImage(varargin{1});
                [varargout{1:nargout}] = residualsPushed(obj, obj.mu, pf, pvox);
            end
            
        case 2
        % ll = obj.residuals(mu, f)
        % mu = varargin{1}
        % f  = varargin{2}
            [varargout{1:nargout}] = residualsWarped(obj, varargin{1}, varargin{2});
            
        case 3
        % ll = obj.residuals(iphi, mu, f, (type))
        % ipsi = varargin{1}
        % mu   = varargin{2}
        % f    = varargin{3}
            if strcmp(type, 'warp')
                wmu = obj.warpTemplate(varargin{1}, varargin{2});
                wmu = obj.reconstructWarpedTemplate(wmu);
                [varargout{1:nargout}] = residualsWarped(obj, wmu, varargin{3});
            else
                [pf, pvox] = obj.pushImage(varargin{3}, varargin{1});
                [varargout{1:nargout}] = residualsPushed(obj, varargin{2}, pf, pvox);
            end
            
        otherwise
            error('Too many arguments');
    end
end
% =========================================================================
function [res, dres] = residualsPushed(obj, mu, f, c)
% FORMAT ([r, dr]) = residualsPushed(obj, mu, f, c)
% obj   - Velocity object
%         Used properties: MatchingTerm, (Normal.s), (Laplace.b)
% mu    - Template image in templatespace
% f     - Observed image pushed in template space
% c     - Pushed voxel count
% (r)   - Residuals [if no argout: write to obj.pres]
% (dr)  - Residuals used for gradient [if no argout: write to obj.presd]
%
% Compute log-likelihood in the template space. Less accurate. Returns:
% 'normal':      res(k)  = 1 / sk * (ck .* muk - fk)^2
%                dres(k) = 1 / sk * (ck .* muk - fk)
% 'laplace':     res(k)  = 1 / bk * abs(ck .* muk - fk)
%                dres(k) = 1 / bk * sign(ck .* muk - fk)
% 'binomial':    res     = f .* log(mu) + (c - f) .* log(1 - mu);
% 'multinomial': res(k)  = sum{k} fk .* log(muk)
    
    % --- Check if nothing to do
    noderiv = {'bernoulli', 'binomial', 'binary', 'categorical', 'multinomial'};
    if nargout == 0 && obj.checkarray('pres') && ...
        (obj.checkarray('presd') || any(strcmpi(obj.MatchingTerm, noderiv)))
        res  = obj.pres;
        dres = obj.presd;
        return
    end
    
    % --- Check if input arrays are usable
    if ~obj.checkarray(mu) || ~obj.checkarray(f) || ~obj.checkarray(c)
        res = [];
        dres = [];
        if obj.Debug
            warning('Cannot compute residuals: missing arrays\n');
        end
        return
    end
    if obj.Debug, fprintf('* residualsPushed\n'); end;
    
    % --- Read dim info
    dim = [size(mu) 1 1 1 1];
    dim = dim(1:4);
    % lat = dim(1:3);
    nc  = dim(4);
    
    % --- Prepare output arrays
    if nargout == 0
        if any(strcmpi(obj.MatchingTerm, {'categorical', 'multinomial'}))
            obj.pres.dim  = dim(1:3);
            obj.presd.dim = dim(1:3);
        else
            obj.pres.dim  = dim;
            obj.presd.dim = dim;
        end
        res  = obj.pres;
        dres = obj.presd;
    else
        if any(strcmpi(obj.MatchingTerm, {'categorical', 'multinomial'}))
            res  = zeros(dim(1:3), 'single');
            dres = zeros(dim(1:3), 'single');
        else
            res  = zeros(dim(1:4), 'single');
            dres = zeros(dim(1:4), 'single');
        end
    end
    
    % --- Deal with different matching terms
    switch lower(obj.MatchingTerm)
        case {'normal', 'gaussian', 'l2'}
            s = obj.Normal.s;   % < sigma^2
            if length(s) == 1
                s = s * ones(1, nc);
            end
            c = single(numeric(c)); % Load count (if on disk, saves i/o)
            for k=1:nc
                % Load one class (if on disk, saves memory)
                mu1 = single(mu(:,:,:,k));
                f1  = single(f(:,:,:,k));
                r = mu1 .* c - f1;
                r(~isfinite(mu1) | ~isfinite(f1) | c == 0) = nan;
                res(:,:,:,k)  = 1/s(k) .* r.^2;
                dres(:,:,:,k) = 1/s(k) .* r;
                clear f1 mu1 r
            end
            
        case {'laplace', 'l1'}
            b = obj.Laplace.b;
            if length(b) == 1
                b = b * ones(1, nc);
            end
            c = single(numeric(c)); % Load count (if on disk, saves i/o)
            for k=1:nc
                % Load one class (if on disk, saves memory)
                mu1 = single(mu(:,:,:,k));
                f1  = single(f(:,:,:,k));
                r = mu1 .* c - f1;
                r(~isfinite(mu1) | ~isfinite(f1) | c == 0) = nan;
                clear mask
                res(:,:,:,k)  = 1/b(k) .* abs(r);
                dres(:,:,:,k) = 1/b(k) .* sign(r);
                clear f1 mu1 r
            end
            
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
            msk = isfinite(mu) & isfinite(f) & c > 0;
            mu(msk) = max(eps('single'), min(1-eps('single'), mu(msk)));
            r = f .* log(mu) + (c - f) .* log(1 - mu);
            r(~msk) = nan;
            res(:,:,:,1) = r;
            clear f mu r
            
        case {'categorical', 'multinomial'}
            % Here, the input template is actually encoded in some kind of  
            % log space (that we will note a). Actual probabilities mu are 
            % mapped back with a softmax function:
            % mu[k] = exp(a[k]) / (sum{l} exp(a[l]))
            
            % E = log p(f | mu, v)
            %   = sum{i,k} [ c(xi) (pf(k)(xi)/c(xi)) log(mu(k)(xi)) ]
            %   = sum{i,k} [ pf(k)(xi) log(mu(k)(xi)) ]
            
            res(:,:,:) = 0;
            for k=1:nc
                % Load one class (if on disk, saves memory)
                mu1 = single(mu(:,:,:,k));
                f1  = single(f(:,:,:,k));
                msk = isfinite(mu) & isfinite(f) & c > 0;
                mu1(msk) = max(eps('single'), min(1-eps('single'), mu1(msk)));
                r = f1 .* mu1;
                r(~msk) = nan;
                res(:,:,:) = res(:,:,:) + r;
                clear f1 mu1 r
            end
        otherwise
            error('Unknown distribution ''%s''', obj.MatchingTerm);
    end
    
    if nargout == 0
        obj.utd.pres = true;
        if any(strcmpi(obj.MatchingTerm, ...
            {'normal', 'gaussian', 'l2', 'laplace', 'l1'}))
            obj.utd.presd = true;
        end
    end
end
% =========================================================================
function res = residualsWarped(obj, mu, f)
% FORMAT (res) = residualsWarped(obj, mu, f)
% obj   - Velocity object
%         Used properties: MatchingTerm, (Normal.s), (Laplace.b)
% mu    - Warped template in image space
% f     - Native image
% (res) - Residuals [if no argout: write to obj.res]
%
% Compute log-likelihood in the native/observed space. More accurate.
% Returns:
% 'normal':      res(k) = 1 / sk * (fk - wmuk)^2
% 'laplace':     res(k) = 1 / bk * abs(fk - wmuk)
% 'binomial':    res    = f .* log(wmu) + (1 - f) .* log(1 - wmu);
% 'multinomial': res(k) = sum{k} fk .* log(wmuk)
        
    % --- Check if nothing to do
    if nargout == 0 && obj.checkarray('res')
        res = obj.res;
        return
    end
    
    % --- Check if input arrays are usable
    if ~obj.checkarray(mu) || ~obj.checkarray(f)
        res = [];
        if obj.Debug
            warning('Cannot compute residuals: missing arrays');
        end
        return
    end
    if obj.Debug, fprintf('* residualsWarped\n'); end;
    
    % --- Read dim info
    dim = [size(mu) 1 1 1 1];
    dim = dim(1:4);
    nc  = dim(4);
    
    % --- Prepare output arrays
    if nargout == 0
        if any(strcmpi(obj.MatchingTerm, {'categorical', 'multinomial'}))
            obj.res.dim  = dim(1:3);
        else
            obj.res.dim  = dim;
        end
        res  = obj.res;
    else
        if any(strcmpi(obj.MatchingTerm, {'categorical', 'multinomial'}))
            res  = zeros(dim(1:3), 'single');
        else
            res  = zeros(dim(1:4), 'single');
        end
    end
    
    % --- Deal with different matching terms
    switch obj.MatchingTerm
        case {'normal', 'gaussian', 'l2'}
            s = obj.Normal.s;   % < sigma^2
            if length(s) == 1
                s = s * ones(1, nc);
            end
            for k=1:nc
                % Load one class (if on disk, saves memory)
                mu1 = single(mu(:,:,:,k));
                f1  = single(f(:,:,:,k));
                r = 1/s(k) .* (f1 - mu1).^2;
                r(~isfinite(mu1) | ~isfinite(f1)) = nan;
                res(:,:,:,k) = r;
                clear f1 mu1 r
            end
            
        case {'laplace', 'l1'}
            b = obj.Laplace.b;
            if length(b) == 1
                b = b * ones(1, nc);
            end
            for k=1:nc
                % Load one class (if on disk, saves memory)
                mu1 = single(mu(:,:,:,k));
                f1  = single(f(:,:,:,k));
                r = 1/b(k) .* abs(f1 - mu1);
                r(~isfinite(mu1) | ~isfinite(f1)) = nan;
                res(:,:,:,k) = r;
                clear f1 mu1 r
            end
            
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
            msk = isfinite(mu) & isfinite(f);
            mu(msk) = max(eps('single'), min(1-eps('single'), mu(msk)));
            r = f .* log(mu) + (1 - f) .* log(1 - mu);
            r(~msk) = nan;
            res(:,:,:,1) = r;
            clear f mu r
            
        case {'categorical', 'multinomial'}
            % Here, the input template is actually encoded in some kind of  
            % log space (that we will note a). Actual probabilities mu are 
            % mapped back with a softmax function:
            % mu[k] = exp(a[k]) / (sum{l} exp(a[l]))
            
            % E = log p(f | mu, v)
            %   = sum{i,k} [ f(k)(xi) log(mu(k)(xi)) ]
            
            res(:) = 0;
            for k=1:nc
                mu1 = single(mu(:,:,:,k));
                f1  = single(f(:,:,:,k));
                msk = isfinite(mu1) & isfinite(f1);
                mu1(msk) = max(eps('single'), min(1-eps('single'), mu1(msk)));
                r = f1 .* log(mu1);
                r(msk) = nan;
                res(:,:,:) = res(:,:,:) + r;
                clear f1 mu1 msk
            end
    end
    
    if nargout == 0
        obj.utd.res = true;
        obj.statusChanged('res');
    end
end