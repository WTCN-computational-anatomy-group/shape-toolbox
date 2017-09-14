function gmu = computeTemplateGrad(obj, mu)
% FORMAT (gmu) = obj.computeTemplateGrad((mu))
% (mu)  - Template [default: obj.mu (normal, laplace) or obj.a (bi, multi)]
% (gmu) - Template spatial gradients [if no argout: write to obj.gmu]
% obj   - Velocity object. This function uses property 'Interpolation'.
%
% Compute spatial gradients (w.r.t. local deformations) of the template
% image.

    % --- Check if nothing to do
    if nargin == 1 && obj.checkarray('gmu')
        gmu = obj.gmu;
        return
    end

    % --- Default arguments
    if nargin < 2
        mu = obj.a; % There's a redirection towards mu for normal/laplace cases
    end
    
    % -- Check that arrays are ready to be used
    if ~obj.checkarray(mu)
        if obj.Debug
            warning('Cannot compute template gradients: missing arrays\n')
        end
        return
    end
    if obj.Debug, fprintf('* computeTemplateGrad\n'); end;

    % --- Read dimensions
    dim = [size(mu) 1 1 1 1];
    dim = dim(1:4);
    lat = dim(1:3);  % Template lattice dimension
    nc  = dim(4);    % Number of classes/modalities
    id  = single(transfo('idmap', lat)); % Identity transform

    % --- Reserve space on disk for the gradient images
    if nargout == 0
        obj.gmu.dim = [dim 3];
        gmu = obj.gmu;
    else
        gmu = zeros([dim 3], 'single');
    end

    % --- Do it
    for i=1:nc
        % Load one class
        mu1 = single(mu(:,:,:,i));
        % Interpolation (Trilinear/Circulant)
        c = spm_diffeo('bsplinc', mu1, obj.Interpolation);
        % Get gradients
        [~, G1, G2, G3] = spm_diffeo('bsplins', c, id, obj.Interpolation);
        % Write data
        gmu(:,:,:,i,1) = G1;
        gmu(:,:,:,i,2) = G2;
        gmu(:,:,:,i,3) = G3;
    end

    if nargout == 0
        obj.utd.gmu = true;
        obj.statusChanged('gmu');
    end
        
end