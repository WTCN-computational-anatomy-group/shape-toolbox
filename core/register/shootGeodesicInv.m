function [iphi, ijac] = shootGeodesicInv(v, vs, prm, itgr, bnd, keep_all)

    if nargin < 6
        keep_all = false;
        if nargin < 5
            bnd = 0; % circulant
            if nargin < 4
                itgr = NaN; % educated guess
            end
        end
    end

    % --- Dimension info
    dim = [size(v) 1 1 1 1];
    dim = dim(1:4);
    lattice_dim = dim(1:3);
    
    N = itgr;
    if ~isfinite(N) || N <= 0
        % Number of time steps from an educated guess about how far to move
        N = double(floor(sqrt(max(max(max( sum(v.^2, 4) ))))) + 1);
    end

    if keep_all
        all_iphi = zeros([lattice_dim 3 N+1], 'single');
        if nargout > 1
            all_ijac = zeros([lattice_dim N+1], 'single');
        else
            all_ijac = [];
        end
    else
        all_iphi = [];
        all_ijac = [];
    end
    
    % Inversion kernel
    F = spm_shoot_greens('kernel', double(lattice_dim), double([vs prm]), bnd);
    % Identity transform
    id = single(spm_warps('identity', lattice_dim));
    % Initial momentum (m_0 = L v_0)
    m = spm_diffeo('vel2mom', single(v), double([vs prm]));

    if ~isempty(all_iphi)
        all_iphi(:,:,:,:,1) = id;
        if ~isempty(all_ijac)
            all_ijac(:,:,:,1) = spm_diffeo('def2det', all_iphi(:,:,:,:,1)-id);
        end
    end
    
    % First iteration
    iphi = id - v/N;
    diphi = iphi;
    
    if ~isempty(all_iphi)
        all_iphi(:,:,:,:,2) = iphi;
        if ~isempty(all_ijac)
            all_ijac(:,:,:,2) = spm_diffeo('def2det', all_iphi(:,:,:,:,2)-id);
        end
    end
    
    % Subsequent iterations
    for t=2:N
        % Jacobian of the step
        % dJ_t = Jac( I - dv_t )
        dJ = spm_diffeo('jacobian', single(diphi));
        m1 = zeros(size(m), 'single');
        % Pointwise matrix multiplication
        % tmp = dJ_t' * m_t
        m1(:,:,:,1) = dJ(:,:,:,1,1).*m(:,:,:,1) + dJ(:,:,:,2,1).*m(:,:,:,2) + dJ(:,:,:,3,1).*m(:,:,:,3);
        m1(:,:,:,2) = dJ(:,:,:,1,2).*m(:,:,:,1) + dJ(:,:,:,2,2).*m(:,:,:,2) + dJ(:,:,:,3,2).*m(:,:,:,3);
        m1(:,:,:,3) = dJ(:,:,:,1,3).*m(:,:,:,1) + dJ(:,:,:,2,3).*m(:,:,:,2) + dJ(:,:,:,3,3).*m(:,:,:,3);
        % Push values to their new locations
        % m_{t+1} = dphi_t( dJ_t' * m_t )
        m = spm_diffeo('pushc', m1, id + v/N);
        % Update velocity
        % v_{t+1} = K m_{t+1}
        v = spm_shoot_greens(m, F, double([vs prm]), bnd);
%         v = spm_diffeo('mom2vel', m, [vs prm 2 2]);
        % Update transform
        diphi = id - v/N;
        try
            iphi(:,:,:,1) = spm_diffeo('bsplins', iphi(:,:,:,1)-id(:,:,:,1), diphi, [1 1 1  1 1 1]) + diphi(:,:,:,1);
            iphi(:,:,:,2) = spm_diffeo('bsplins', iphi(:,:,:,2)-id(:,:,:,2), diphi, [1 1 1  1 1 1]) + diphi(:,:,:,2);
            iphi(:,:,:,3) = spm_diffeo('bsplins', iphi(:,:,:,3)-id(:,:,:,3), diphi, [1 1 1  1 1 1]) + diphi(:,:,:,3);
        catch
            error('Unexpected stuff happening')
        end
        
        % save
        if ~isempty(all_iphi)
            all_iphi(:,:,:,:,t+1) = iphi;
            if ~isempty(all_ijac)
                all_ijac(:,:,:,t+1) = spm_diffeo('def2det', all_iphi(:,:,:,:,t+1)-id);
            end
        end
    end
    
    if ~isempty(all_iphi)
        iphi = all_iphi;
        ijac = all_ijac;
    else
        if nargout > 1
            ijac = spm_diffeo('def2det', single(iphi-id));
        end
    end
end