function varargout = exponentiateVelocity(obj, varargin)
% FORMAT ([iphi, iJ, phi, J]) = obj.exponentiateVelocity((v), ('iphi'), ('ijac'), ('phi'), ('jac'))
% Exponentiate the velocity to recover the inverse transforms

    % --- Check if velocity provided
    if numel(varargin) > 0 && ~ischar(varargin{1})
        v = varargin{1};
        varargin = varargin(2:end);
    else
        v = [];
    end
    
    % --- Read which computations to perform
    do      = struct;
    do.phi  = false;
    do.iphi = false;
    do.jac  = false;
    do.ijac = false;
    for i=1:length(varargin)
        do.(lower(varargin{i})) = true;
    end
    if do.ijac
        do.iphi = true;
    end
    if do.jac
        do.phi = true;
    end
    
    % --- Check nothing to do
    out = struct;
    if isempty(v)
        do_something = false;
        for i=1:length(varargin)
            if obj.utd.(lower(varargin{i}))
                out.(lower(varargin{i})) = obj.(lower(varargin{i}));
                do.(lower(varargin{i})) = false;
            end
            do_something = do_something | do.(lower(varargin{i}));
        end
        if ~do_something
            varargout = cell(size(varargin));
            for i=1:length(varargin)
                varargout{i} = out.(lower(varargin{i}));
            end
            return
        end
    end
        
    
    % --- Default parameters
    if isempty(v)
        obj.reconstructVelocity();
        v = obj.v;
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
    if obj.Debug
        outputs = fieldnames(do);
        strout = '';
        for i=1:length(outputs)
            if do.(outputs{i})
                strout = [strout outputs{i} ','];
            end
        end
        strout = strout(1:end-1);
        fprintf('* exponentiateVelocity: %s\n', strout);
    end;
    
    % --- Load data in memory
    if isa(v, 'file_array')
        v = single(numeric(v));
    end
    
    % --- Deal with different cases
    if do.iphi && ~do.phi
        if do.ijac
            [out.iphi, out.ijac] = pushIPhiJac(obj, v);
        else
            out.iphi = pushIPhi(obj, v);
        end
    elseif do.phi && ~do.iphi
        if do.jac
            [out.phi, out.jac] = pushPhiJac(obj, v);
        else
            out.phi = pushPhi(obj, v);
        end
    elseif do.phi && do.iphi
        if (do.jac || do.ijac)
            [out.iphi, out.ijac, out.phi, out.jac] = pushIPhiPhiJac(obj, v);
        else
            [out.iphi, out.phi] = pushIPhiPhi(obj, v);
        end
    end

    % --- Write on disk
    if nargout == 0
        if do.iphi
            obj.iphi.dim = size(out.iphi);
            obj.iphi(:)  = out.iphi(:);
            obj.statusChanged('iphi');
        end
        if do.phi
            obj.phi.dim = size(out.phi);
            obj.phi(:)  = out.phi(:);
            obj.statusChanged('phi');
        end
        if do.ijac
            obj.ijac.dim = size(out.ijac);
            obj.ijac(:)  = out.ijac(:);
            obj.statusChanged('ijac');
        end
        if do.jac
            obj.jac.dim = size(out.jac);
            obj.jac(:)  = out.jac(:);
            obj.statusChanged('jac');
        end
    else
        varargout = cell(size(varargin));
        for i=1:length(varargin)
            varargout{i} = out.(lower(varargin{i}));
        end
    end
end

% === Specialized/Optimized exonentiators =================================

function phi = pushPhi(obj, v)
    phi = spm_shoot3d(v, [obj.VoxelSize obj.RegParam], obj.Integration);
end

function [phi, jac] = pushPhiJac(obj, v)
    [phi, jac] = spm_shoot3d(v, [obj.VoxelSize obj.RegParam], obj.Integration);
end

% function iphi = pushIPhi(obj, v)
%     iphi = spm_shoot3di(v, [obj.VoxelSize obj.RegParam], obj.Integration);
% end

function [iphi, ijac] = pushIPhiJac(obj, v)
    [iphi, ijac] = spm_shoot3di(v, [obj.VoxelSize obj.RegParam], obj.Integration);
end

function [iphi, phi] = pushIPhiPhi(obj, v)
    [phi, ~, ~, iphi, ~] = spm_shoot3d(v, [obj.VoxelSize obj.RegParam], obj.Integration);
end

function [iphi, ijac, phi, jac] = pushIPhiPhiJac(obj, v)
    [phi, jac, ~, iphi, ijac] = spm_shoot3d(v, [obj.VoxelSize obj.RegParam], obj.Integration);
end

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