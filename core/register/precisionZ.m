function ww = precisionZ(w, varargin)
% FORMAT ww = precisionZ(w, (vs), (prm))
% w   - Principal subspace
% vs  - Voxel size of the initial velocity lattice [1 1 1]
% prm - Parameters of the L operator (see spm_diffeo)
%       [0.0001 0.001 0.2 0.05 0.2]
% ww  - Precision matrix of P(Z | W), i.e., W' * L * W
%
% Compute the precision matrix of p(z|W), where Z are the latent
% coordinates of the velocity.
% (Also named 'inverse covariance' or 'regularization' matrix)
    
    % --- Parse inputs
    p = inputParser;
    p.FunctionName = 'precisionZ';
    p.addRequired('w', @checkarray);
    p.addOptional('vs', [1 1 1]);
    p.addOptional('prm', [0.0001 0.001 0.2 0.05 0.2]);
    p.addParameter('output', []);
    p.addParameter('debug', false);
    p.parse(w, varargin{:});
    
    if p.Results.debug, fprintf('* precisionZ\n'); end;
    
    % --- Dim info
    dim         = [size(w) 1 1 1];
    dim         = dim(1:5);
    dim_latent  = dim(5);       % Number of principal components

    % --- Init output
    % To save i/o, we do computations in memory and write on disk at the
    % end
    ww = zeros(dim_latent, 'double');
    
    % --- Compute WW
    % 1) Compute momentum of the k1-th basis vector (L * Wk1)
    % 2) Left-multiply by the k2-the basis vector (Wk2' * L * Wk1)
    for k1=1:dim_latent
        w1 = numeric(w(:,:,:,:,k1));
        u = spm_diffeo('vel2mom', single(w1), [p.Results.vs p.Results.prm]);
        ww(k1,k1) = sumall(u .* double(w1));
        clear w1
        for k2=(k1+1):dim_latent
            ww(k1, k2) = sumall(u .* double(w(:,:,:,:,k2)));
            ww(k2, k1) = ww(k1, k2);
        end
    end
        
    % --- Write on disk
    if ~isempty(p.Results.output)
        ww = saveOnDisk(p.Results.output, ww, 'name', 'ww');
    end
        
end