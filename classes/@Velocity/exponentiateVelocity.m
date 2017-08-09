function varargout = exponentiateVelocity(obj, varargin)
% FORMAT ([iphi, iJ, phi, J]) = obj.exponentiateVelocity((v), ('iphi'), ('ijac'), ('phi'), ('jac'))
% Exponentiate the velocity to recover the inverse transforms

    % --- Check if velocity provided
    if numel(varargin) > 0 && ~ischar(varargin{1})
        v = varargin{1};
        varargin = varargin(2:end);
    else
        if ~obj.checkarray('v')
            obj.reconstructVelocity();
        end
        v = obj.v;
    end
    
    % --- Read which computations to perform
    do      = struct;
    do.phi  = false;
    do.iphi = true;
    do.jac  = false;
    do.ijac = false;
    for i=1:length(varargin)
        do.(varargin{i}) = true;
    end
    
    % --- If input not usable, exit
    if ~obj.checkarray(v)
        varargout{1} = single([]);
        varargout{2} = single([]);
        varargout{3} = single([]);
        varargout{4} = single([]);
        if obj.Debug
            warning('Cannot exponentiate velocity: missing arrays');
        end
        return;
    end
    if obj.Debug, fprintf('* exponentiateVelocity\n'); end;
    
    % --- Load data in memory
    if isa(v, 'file_array')
        v = single(numeric(v));
    end
    
    % --- Deal with different cases
    if do.iphi
        if ~do.ijac && ~do.phi && ~do.jac
            if nargout == 0 && obj.utd.iphi,  return; end % Nothing to do
            iphi = pushIPhi(obj, v);
            varargout{1} = iphi;
        end
    end

    % --- Write on disk
    if nargout == 0
        if do.iphi
            obj.iphi.dim = size(iphi);
            obj.iphi(:)  = iphi(:);
            obj.utd.iphi = true;
            obj.statusChanged('iphi');
        end
    end
    
    % For now I only do iphi because that's everything we need for the
    % registration. What would probably be better would be to create a
    % 'classical' Velocity class which implements geodesic shooting, and
    % embed it in this 'latent' Velocity class.
    % At the same time, if W = Id, the latent and clasic classes are
    % identical, so it makes sense to leave everything at the same place.
end

% === Specialized/Optimized exonentiators =================================

function iphi = pushIPhi(obj, v)
% FORMAT iphi = pushIPhi(obj, v, (N), (odim))
% v      - Initial velocity (in template space)
% iphi   - Inverse transform (warps template to image)
%
% If only iphi is needed, exponentiate by pushing.

    % --- Dimension info
    dim = [size(v) 1 1 1 1];
    dim = dim(1:4);
    lattice_dim = dim(1:3);
    
    N = obj.Integration;
    if ~isfinite(N) || N <= 0
        % Number of time steps from an educated guess about how far to move
        N = double(floor(sqrt(max(max(max( sum(v.^2, 4) ))))) + 1);
    end

    % Inversion kernel
    F = spm_shoot_greens('kernel', lattice_dim, [obj.VoxelSize obj.RegParam]);
    % Identity transform
    id = single(transfo('idmap', lattice_dim));
    % Initial momentum (m_0 = L v_0)
    m = spm_diffeo('vel2mom', v, [obj.VoxelSize obj.RegParam]);

    % First iteration
    iphi = id - v/N;
    diphi = iphi;
    if obj.Verbose, fprintf('Shoot |>'); end;
    
    % Subsequent iterations
    for t=2:N
        % Jacobian of the step
        % dJ_t = Jac( I - dv_t )
        dJ = spm_diffeo('jacobian', diphi);
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
        v = spm_shoot_greens(m, F, [obj.VoxelSize obj.RegParam]);
        % Update transform
        diphi = id - v/N;
        iphi(:,:,:,1) = spm_diffeo('bsplins', iphi(:,:,:,1)-id(:,:,:,1), diphi, [1 1 1  1 1 1]) + diphi(:,:,:,1);
        iphi(:,:,:,2) = spm_diffeo('bsplins', iphi(:,:,:,2)-id(:,:,:,2), diphi, [1 1 1  1 1 1]) + diphi(:,:,:,2);
        iphi(:,:,:,3) = spm_diffeo('bsplins', iphi(:,:,:,3)-id(:,:,:,3), diphi, [1 1 1  1 1 1]) + diphi(:,:,:,3);
        % Verbose
        if obj.Verbose, fprintf('>'); end;
    end
    if obj.Verbose, fprintf('| %d integration steps\n', N); end;
end