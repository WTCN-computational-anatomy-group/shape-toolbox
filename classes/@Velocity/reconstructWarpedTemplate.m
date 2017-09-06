function wmu = reconstructWarpedTemplate(obj, wa)
% FORMAT (wmu) = obj.reconstructTemplate((wa))
% (wa)  - Warped log-Template [default: obj.wa]
% (wmu) - Warped template [if no argout: write to obj.wmu]
% obj   - Velocity object. This function uses property 'MatchingTerm'.
%
% Reconstruct the template probability maps from their log-space. Only
% useful for matching terms binomial and multinomial.

    % --- Check if nothing to do
    if nargout == 0 && obj.checkarray('wmu')
        return
    end
    
    % --- Default arguments
    if nargin < 2
        obj.warpTemplate();
        wa = obj.wa;
    end
    if ~any(strcmpi(obj.MatchingTerm, ...
            {'bernoulli', 'binomial', 'binary', 'categorical', 'multinomial'}))
        wmu = wa;
        return
    end
    
    % -- Check that arrays are ready to be used
    if ~obj.checkarray(wa)
        if obj.Debug
            warning('Cannot reconstruct template: missing arrays\n')
        end
        return
    end
    if obj.Debug, fprintf('* reconstructTemplate\n'); end;

    % --- Read dimensions
    dim = [size(wa) 1 1 1 1];
    dim = dim(1:4);
    lat = dim(1:3);  % Template lattice dimension
    nc  = dim(4);    % Number of classes/modalities

    % --- Reserve space on disk for the gradient images
    if nargout == 0
        obj.wmu.dim = dim;
        wmu = obj.wmu;
    else
        wmu = zeros(dim, 'single');
    end

    % --- Do it
    switch lower(obj.MatchingTerm)
        case {'bernoulli', 'binomial', 'binary'}
            wa = exp(single(numeric(wa)));
            wmu(:,:,:) = wa ./ (1 + wa);
        case {'categorical', 'multinomial'}
            s = zeros(lat, 'single');
            for k=1:nc
                ea = exp(single(wa(:,:,:,k)));
                wmu(:,:,:,k) = ea;
                s = s + ea;
            end
            for k=1:nc
                wmu(:,:,:,k) = wmu(:,:,:,k) ./ s;
            end
    end

    if nargout == 0
        obj.statusChanged('wmu');
        obj.utd.wmu = true;
    end
        
end