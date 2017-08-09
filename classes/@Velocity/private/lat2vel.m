function v = lat2vel(z, w, v)
% FORMAT v = lat2vel(z, w, (v))
% z   - (array or file_array) Latent coordinates [K (n_obs)]
% w   - (array or file_array) Basis functions    [nx ny nz 3 K]
% (v) - (file_array) Specify disk output (allows to save memory)

    % --- Read dimensions
    dim  = [size(w) 1 1 1 1 1];
    dim = dim(1:5);
    lat = dim(1:3);    % Lattice dimension ([nx ny nz])
    vec = dim(4);      % Vector dimension (should be 3)
    K   = dim(5);      % Latent (low) dimension
    N   = size(z, 2);  % Number of observations (i.e. latent velocities)
    
    % --- Allocate output
    if nargin < 3
        v = zeros([lat vec N], 'single');
    else
        v.dim = [lat vec N]; % Just in case
    end

    % --- Compute linear combination
    if isa(w, 'file_array')
        % Save memory w.r.t. basis functions and observations 
        for n=1:N
            for k=1:K
                v(:,:,:,:,n) = v(:,:,:,:,n) + w(:,:,:,:,k) * z(k,n);
            end
        end
    elseif isa(z, 'file_array')
        % Save memory w.r.t. observations
        for n=1:N
            for k=1:K
                v(:,:,:,:,n) = v(:,:,:,:,n) + w * z(:,n);
            end
        end
    else
        % Save time
        v(:,:,:,:,:) = v(:,:,:,:,:) + reshape(w, [], latent_dim) * z;
    end
    
end