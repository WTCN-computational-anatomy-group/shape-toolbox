function t = latentToTensorCov(w, h)

    % --- Dim info
    dim         = [size(w) 1 1 1];
    dim         = dim(1:5);
    dim_lattice = dim(1:3);
    dim_vec     = dim(4);
    dim_latent  = dim(5);
    
    % --- Allocate output
    t = zeros([dim_lattice dim_vec dim_vec], 'single');
    h = numeric(h);
    
    % --- Do it
    for v1=1:dim_vec
        for v2=1:dim_vec
            acc = zeros(dim_lattice, 'double');
            for k=1:dim_latent
                wk = w(:,:,:,v1,k);
                for l=1:dim_latent
                    wl = w(:,:,:,v2,l);
                    acc = acc + wk .* h(k,l) .* wl;
                end
            end
            t(:,:,:,v1,v2) = acc;
        end
    end
            

end