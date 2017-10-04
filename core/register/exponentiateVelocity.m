function varargout = exponentiateVelocity(v, varargin)
% FORMAT ([iphi, iJ, phi, J]) = exponentiateVelocity(v,
%   ('iphi'), ('ijac'), ('phi'), ('jac'),
%   ('itgr', itgr),
%   ('vs',   vs),
%   ('prm',  prm),
%
% ** Required **
% v      - Inital velocity (in template space)
% ** Optional **
% 'iphi' - Compute iphi (f = mu(iphi) -> Warps template to image)
% 'ijac' - Compute D(iphi)
% 'phi'  - Compute phi (f(phi) = mu -> Warps image to template)
% 'jac'  - Compute D(phi)
% ** Keyword arguments ии
% itgr   - Number of integration step (for geodesic shooting) [auto]
% vs     - Voxel size of the initial velocity lattice [1 1 1]
% prm    - Parameters of the L operator (see spm_diffeo) 
%          [0.0001 0.001 0.2 0.05 0.2]
%
% Exponentiate the initial velocity to recover transforms.

    % --- Read which computations to perform
    do      = struct;
    do.phi  = false;
    do.iphi = false;
    do.jac  = false;
    do.ijac = false;
    which   = {};
    i = 1;
    while ~isempty(varargin) && i <= numel(varargin)
        if ischar(varargin{i}) && any(strcmpi(varargin{i}, {'iphi', 'phi', 'ijac', 'jac'}))
            which{end+1} = lower(varargin{i});
            do.(lower(varargin{i})) = true;
            varargin = {varargin{1:i-1} varargin{i+1:end}};
        else
            i = i + 1;
        end
    end
    if do.ijac
        do.iphi = true;
    end
    if do.jac
        do.phi = true;
    end
    
    % --- If v is empty -> return empty
    if isempty(v)
        varargout = cell(size(which));
        for i=1:numel(which)
            varargout{i} = [];
        end
        return
    end
    
    % --- Parse inputs
    p = inputParser;
    p.FunctionName = 'exponentiateVelocity';
    p.addRequired('v');
    p.addParameter('itgr',   nan);
    p.addParameter('vs',     [1 1 1]);
    p.addParameter('prm',    [0.0001 0.001 0.2 0.05 0.2]);
    p.addParameter('output', []);
    p.addParameter('debug',  false);
    p.parse(v, varargin{:});
    itgr   = p.Results.itgr;
    vs     = p.Results.vs;
    prm    = p.Results.prm;
    debug  = p.Results.debug;
    output = p.Results.output;
    
    
    if debug, fprintf('* exponentiateVelocity\n'); end;
    
    % --- Load data in memory
    v = single(numeric(v));
    
    % --- Deal with different cases
    if do.iphi && ~do.phi
        if do.ijac
            [out.iphi, out.ijac] = pushIPhiJac(v, itgr, vs, prm);
        else
            out.iphi = pushIPhi(v, itgr, vs, prm);
        end
    elseif do.phi && ~do.iphi
        if do.jac
            [out.phi, out.jac] = pushPhiJac(v, itgr, vs, prm);
        else
            out.phi = pushPhi(v, itgr, vs, prm);
        end
    elseif do.phi && do.iphi
        if (do.jac || do.ijac)
            [out.iphi, out.ijac, out.phi, out.jac] = pushIPhiPhiJac(v, itgr, vs, prm);
        else
            [out.iphi, out.phi] = pushIPhiPhi(v, itgr, vs, prm);
        end
    end

    % --- Write on disk // Return objects in the right order
    varargout = cell(size(which));
    if ~iscell(output)
        output = {output};
    end
    for i=1:numel(which)
        if numel(output) >= i && ~isempty(output{i})
            varargout{i} = saveOnDisk(output{i}, out.(which{i}), 'name', which{i});
        else
            varargout{i} = out.(which{i});
        end
    end
end

% === Specialized/Optimized exponentiators ================================

function phi = pushPhi(v, itgr, vs, prm)
    phi = spm_shoot3d(v, [vs prm], itgr);
end

function [phi, jac] = pushPhiJac(v, itgr, vs, prm)
    [phi, jac] = spm_shoot3d(v, [vs prm], itgr);
end

function [iphi, ijac] = pushIPhiJac(v, itgr, vs, prm)
    [iphi, ijac] = pushIPhi(v, itgr, vs, prm);
end

function [iphi, phi] = pushIPhiPhi(v, itgr, vs, prm)
    [phi, ~, ~, iphi, ~] = spm_shoot3d(v, [vs prm], itgr);
end

function [iphi, ijac, phi, jac] = pushIPhiPhiJac(v, itgr, vs, prm)
    [phi, jac, ~, iphi, ijac] = spm_shoot3d(v, [vs prm], itgr);
end

function [iphi, ijac] = pushIPhi(v, itgr, vs, prm)

    % --- Dimension info
    dim = [size(v) 1 1 1 1];
    dim = dim(1:4);
    lattice_dim = dim(1:3);
    
    N = itgr;
    if ~isfinite(N) || N <= 0
        % Number of time steps from an educated guess about how far to move
        N = double(floor(sqrt(max(max(max( sum(v.^2, 4) ))))) + 1);
    end

    % Inversion kernel
    F = spm_shoot_greens('kernel', lattice_dim, [vs prm]);
    % Identity transform
    id = single(transfo('idmap', lattice_dim));
    % Initial momentum (m_0 = L v_0)
    m = spm_diffeo('vel2mom', v, [vs prm]);

    % First iteration
    iphi = id - v/N;
    diphi = iphi;
    
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
        v = spm_shoot_greens(m, F, [vs prm]);
        % Update transform
        diphi = id - v/N;
        iphi(:,:,:,1) = spm_diffeo('bsplins', iphi(:,:,:,1)-id(:,:,:,1), diphi, [1 1 1  1 1 1]) + diphi(:,:,:,1);
        iphi(:,:,:,2) = spm_diffeo('bsplins', iphi(:,:,:,2)-id(:,:,:,2), diphi, [1 1 1  1 1 1]) + diphi(:,:,:,2);
        iphi(:,:,:,3) = spm_diffeo('bsplins', iphi(:,:,:,3)-id(:,:,:,3), diphi, [1 1 1  1 1 1]) + diphi(:,:,:,3);
    end
    
    if nargout > 1
        ijac = spm_diffeo('jacobian', iphi);
    end
end