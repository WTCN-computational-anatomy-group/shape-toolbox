function ww = computeRegZ(obj, w)
% FORMAT (regz) = obj.computeRegZ((w))
% (w)  - Principal subspace [default: obj.w]
% (ww) - Inverse covariance matrix of P(Z | W), i.e., W' * L * W
% obj  - Velocity object. This function uses properties 'VoxelSize' and
%        'RegParam'
%
% Compute the inverse covariance matrix of p(Z|W) (regularization on Z),
% i.e., W' * L * W.


    % --- Check up-to-date
    if nargin == 1 && obj.checkarray('regz')
        ww = obj.regz;
        return
    end

    % --- Default arguments
    if nargin < 2
        w = obj.w;
    end
    
    % --- Check that all arrays are ready to be used
    if ~obj.checkarray(w)
        ww = [];
        if obj.Debug
            warning('Cannot compute Z regularization matrix. Missing arrays');
        end
        return
    end
    if obj.Debug, fprintf('* computeRegZ\n'); end;
    
    % --- Dim info
    dim         = [size(w) 1 1 1];
    dim         = dim(1:5);
%     dim_lattice = dim(1:3);     % Dimension of the initial velocity lattice
%     dim_vec     = dim(4);       % Dimension of the field (2/3)
    dim_latent  = dim(5);       % Number of principal components

    % --- Init output
    % To save i/o, we do computations in memory and write on disk at the
    % end
    ww = zeros(dim_latent, 'double');
    
    % --- Compute WW
    % The way it is currently done, if w is on disk, there is
    % (dim_latent^2) i/o operations. Maybe we could save them, bu then we
    % need to load the full W on memory. Unless we split along Z-slices ?
    %
    % However, it's not something which should (usually) be done here, it 
    % is only done outside of the velocity-loop by the PrincipalGeodesic
    % object.
    for k1=1:dim_latent
        % Compute momentum of the k1-th basis vector (L * Wk1)
        u = spm_diffeo('vel2mom', w(:,:,:,:,k1), [obj.VoxelSize obj.RegParam]);
        % Left-multiply by the k2-the basis vector (Wk2' * L * Wk1)
        for k2=k1:dim_latent
            ww(k1, k2) = sumall(u .* w(:,:,:,:,k2));
            ww(k2, k1) = ww(k1, k2);
        end
    end
    
    % --- Write on disk
    if nargout == 0
        obj.regz.dim    = size(ww);
        obj.regz(:)     = ww(:);
        ww              = obj.regz;
        obj.statusChanged('regz');
    end
        
end