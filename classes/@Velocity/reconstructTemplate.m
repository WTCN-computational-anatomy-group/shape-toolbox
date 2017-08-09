function mu = reconstructTemplate(obj, a)
% FORMAT (mu) = obj.reconstructTemplate((a))
% (a)   - Log-Template [default: obj.a]
% (mu)  - Template [if no argout: write to obj.mu]
% obj   - Velocity object. This function uses property 'MatchingTerm'.
%
% Reconstruct the template probability maps from their log-space. Only
% useful for matching terms binomial and multinomial.

    % --- Check if nothing to do
    if nargout == 0 && obj.checkarray('mu')
        return
    end
    if ~any(strcmpi(obj.MatchingTerm, ...
            {'bernoulli', 'binomial', 'binary', 'categorical', 'multinomial'}))
        return
    end

    % --- Default arguments
    if nargin < 2
        a = obj.a;
    end
    
    % -- Check that arrays are ready to be used
    if ~obj.checkarray(a)
        if obj.Debug
            warning('Cannot reconstruct template: missing arrays\n')
        end
        return
    end
    if obj.Debug, fprintf('* reconstructTemplate\n'); end;

    % --- Read dimensions
    dim = [size(a) 1 1 1 1];
    dim = dim(1:4);
    lat = dim(1:3);  % Template lattice dimension
    nc  = dim(4);    % Number of classes/modalities

    % --- Reserve space on disk for the gradient images
    if nargout == 0
        obj.mu.dim = dim;
        mu = obj.mu;
    else
        mu = zeros(dim, 'single');
    end

    % --- Do it
    switch lower(obj.MatchingTerm)
        case {'bernoulli', 'binomial', 'binary'}
            a = exp(single(numeric(a)));
            mu(:,:,:) = a ./ (1 + a);
        case {'categorical', 'multinomial'}
            s = zeros(lat, 'single');
            for k=1:nc
                ea = exp(single(a(:,:,:,k)));
                mu(:,:,:,k) = ea;
                s = s + ea;
            end
            for k=1:nc
                mu(:,:,:,k) = mu(:,:,:,k) ./ s;
            end
    end

    if nargout == 0
        obj.statusChanged('mu');
        obj.utd.mu = true;
    end
        
end