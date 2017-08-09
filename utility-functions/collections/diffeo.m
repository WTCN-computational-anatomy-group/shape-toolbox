function varargout = diffeo(id, varargin)
% FORMAT L = diffeo('penalty', lat_dim, lat_vs, param)
% FORMAT K = diffeo('kernel', L, lat_dim)
% FORMAT K = diffeo('greens', L, lat_dim)
% FORMAT m = diffeo('vel2mom', v, L, lat_dim)
% FORMAT v = diffeo('mom2vel', m, K, lat_dim)
% FORMAT [iT, ...] = diffeo('pullexp', v, vs, N, L, K)
% FORMAT [T, ...]  = diffeo('pushexp',v, vs, N, L, K)
%
% Collection of tools to compute diffeomorphic transformations.
%
% FORMAT help diffeo>function
% Returns the help file of the selected function.
    switch id
        case 'penalty'
            [varargout{1:nargout}] = penalty(varargin{:});
        case 'kernel'
            [varargout{1:nargout}] = kernel(varargin{:});
        case 'greens'
            [varargout{1:nargout}] = greens(varargin{:});
        case 'vel2mom'
            [varargout{1:nargout}] = vel2mom(varargin{:});
        case 'mom2vel'
            [varargout{1:nargout}] = mom2vel(varargin{:});
        case 'pullexp'
            [varargout{1:nargout}] = pullexp(varargin{:});
        case 'pushexp'
            [varargout{1:nargout}] = pushexp(varargin{:});
        otherwise
            error('No function named %s\n', string(id));
    end
end

%% Functions: PULLEXP/PUSHEXP

function [T, J] = pushexp(v, vs, N, L, K)
% FORMAT [T, J] = diffeo('pushexp',v, vs, N, L, (K))
% FORMAT [T, J] = diffeo('pushexp',v, vs, N, param)
% v     - Initial velocity field
% vs    - Lattice voxel size
% N     - Number of integration steps
% L     - Regularization matrix obtained with diffeo('penalty')
% K     - Greens function obtained with diffeo('greens')
% param - Regularization parameters. They can be provided instead of L and
%         K. These matrices will have to be computed anyway, so this format
%         might be slower.
% T     - Direct transform (T(source) = target)
% J     - Jacobian of the direct transform
%
% Exponentiate the initial velocity field (= compute the corresponding 
% transformation) by "pushing" the forward transform.
% Equations (17-20) of Ashburner's geodesic shooting paper.
% See: Ashburner and Friston, "Diffeomorphic registration using geodesic
% shooting aud Gauss-Newton optimisation", NeuroImage (2011).

    % --- Parameters ---
    dim = size(v);              % Dimensions of the velocity volume
    lat_dim = dim(1:(end-1));   % Dimensions of the lattice
    vec_dim = dim(end);         % Should be ~ im_dim = length(lat_dim)
    if nargin < 5 && issame(size(L), [1 5])
        % Parameters are provided, not matrices
        param = L;
        L = penalty(lat_dim, vs, param);
        K = greens(L, lat_dim);
    elseif nargin < 5
        % For some reason, L is provided but not K
        K = greens(L, lat_dim);
    elseif issame(size(L), [1 5])
        % For some reason, K is provided but not L
        param = L;
        L = penalty(lat_dim, vs, param);
    end
    
    % --- Useful constants ---
    % Store identity map
    I = transfo('idmap', lat_dim, vs);
    % Store jacobian of identity
    one_jac = diag(1./vs);
    one_jac = reshape(one_jac, [ones([1 vec_dim]) vec_dim vec_dim]);
    JI = repmat(one_jac, lat_dim);
    clear one_jac
    % Store jacobian operator
    fJ = moperators('jacobian', lat_dim, vs, vec_dim);
    
    % --- Initialize ---
    % Compute m0 (eq.8)
    m0 = diffeo('vel2mom', v, L);
    % Initialise the forward transforms at identity (eq.10)
    T = I;
    % Initialise the Jacobians at identity (eq.10)
    J = JI;
    
    % --- Euler-like integration ---
    fprintf('start |');
    for n = 1:N
        fprintf('-');
                
        % Update transformation (eq.17)
        % T_{t+dt} = dT o T_{t}
        % with dT = I + v_{t}/N
        T = transfo('compose', I + v/N, vs, T, vs);
        
        % Update jacobians (eq.18)
        lh = fJ * arrayutils('field2v', v/N);   % <- Jac(dv)
        lh = reshape(lh, [lat_dim vec_dim vec_dim]);
        lh = pointwise_expm(lh);                % <- Exp(Jac(dv))
        lh = pull_jacobian(lh, T, vs);          % <- Exp(Jac(dv)) o T_{t}
        J = pointwise_multrans_jacobian(lh, J); % <- (Exp(Jac(dv)) o T_{t}) '* J_{t}
        clear lh
        
        % Update momentum (eq.19)
        iJ = pointwise_invert(J);                % <- J^{-1}
        m = pointwise_multrans_jacobian(iJ, m0); % <- J^{-1}' * u0
        clear iJ
        m = push_field(T, m, vs);                % <- T*(J^{-1}' * u0)
        
        % Update velocity (eq.16)
        % v_{t+dt} = K u_{t+dt}
        v = diffeo('mom2vel', m, K);
    end
    fprintf('| end\n');
end

function [iT, iJ, T, J] = pullexp(v, vs, N, L, K)
% FORMAT [iT, iJ, (T, J)] = diffeo('pullexp', v, vs, N, L, (K))
% FORMAT [iT, iJ, (T, J)] = diffeo('pullexp',v, vs, N, param)
% v     - Initial velocity field
% vs    - Lattice voxel size
% N     - Number of integration steps
% L     - Regularization matrix obtained with diffeo('penalty')
% K     - Greens function obtained with diffeo('greens')
% param - Regularization parameters. They can be provided instead of L and
%         K. These matrices will have to be computed anyway, so this format
%         might be slower.
% iT    - Inverse transform (iT(target) = source)
% iJ    - Jacobian of the inverse transform
% (T)   - Direct transform (T(source) = target)
% (J)   - Jacobian of the direct transform
%
% Exponentiate the initial velocity field (= compute the corresponding 
% transformation) by "pulling" the inverse transform.
% Equations (8-16) of Ashburner's geodesic shooting paper.
% See: Ashburner and Friston, "Diffeomorphic registration using geodesic
% shooting aud Gauss-Newton optimisation", NeuroImage (2011).
    
    % --- Parameters ---
    dim = size(v);
    lat_dim = dim(1:(end-1));
    vec_dim = dim(end);
    if nargin < 5 && issame(size(L), [1 5])
        % Parameters are provided, not matrices
        param = L;
        L = penalty(lat_dim, vs, param);
        K = greens(L, lat_dim);
    elseif nargin < 5
        % For some reason, L is provided but not K
        K = greens(L, lat_dim);
    elseif issame(size(L), [1 5])
        % For some reason, K is provided but not L
        param = L;
        L = penalty(lat_dim, vs, param);
    end
    
    % --- Useful constants ---
    % Store identity map
    I = transfo('idmap', lat_dim, vs);
    % Store jacobian of identity
    one_jac = diag(1./vs);
    one_jac = reshape(one_jac, [ones([1 vec_dim]) vec_dim vec_dim]);
    JI = repmat(one_jac, lat_dim);
    clear one_jac
    % Store jacobian operator
    fJ = moperators('jacobian', lat_dim, vs, vec_dim);
    
    % --- Initialize ---
    % Compute m0 (eq.8)
    m0 = diffeo('vel2mom', v, L);
    % Initialise the direct and inverse transforms at identity (eq.9-10)
    iT = I;
    if nargout > 2
        T = I;
    end
    % Initialise the Jacobians at identity (eq.9-10)
    iJ = JI;
    if nargout > 2
        J = JI;
    end
    
    % --- Euler-like integration ---
    fprintf('start |');
    for n = 1:N
        fprintf('-');
        
        % Update inverse transformation (eq.11)
        % iT_{t+dt} = iT_{t} o diT
        % with diT = I - v_{t}/N
        iT = transfo('compose', iT, vs, I - v/N, vs);
        
        % Update inverse jacobian (eq.12)
        % iJ_{t+dt} = (iJ_{t} o diT)' * Jac(diT)
        % (o and '* are pointwise operations, i.e., matrix 
        % multiplications at each point of the lattice)
        iJ = pull_jacobian(iJ, I - v/N, vs);              % <- iJ_{t} o diT
        rh = fJ * arrayutils('field2v', v/N);             % <- Jac(dv)
        rh = JI - reshape(rh, [lat_dim vec_dim vec_dim]); % <- Jac(diT) = Jac(I) - Jac(dv)
        iJ = pointwise_multrans_jacobian(iJ, rh);         % <- iJ_{t+dt}
        clear rh;
        
        if nargout > 2
            % Update direct transformation (eq.13)
            % T_{t+dt} = dT o T_{t}
            % with dT = I + v_{t}/N
            T = transfo('compose', I + v/N, vs, T, vs);
            
            % Update direct jacobian (eq.14)
            % J_{t+dt} = (Jac(dT) o T_{t})' * J_{t}
            % (o and '* are pointwise operations, i.e., matrix 
            % multiplications at each point of the lattice)
            lh = fJ * arrayutils('field2v', v/N);             % <- Jac(dv)
            lh = JI + reshape(lh, [lat_dim vec_dim vec_dim]); % Jac(dT) = Jac(I) + Jac(dv)
            lh = pull_jacobian(lh, T, vs);                    % <- Jac(dT) o T_{t}
            J = pointwise_multrans_jacobian(lh, J);           % <- J_{t+dt}
            clear lh
        end
        
        % Update momentum (eq.15)
        % u_{t+dt} = det(iJ_{t-dt}) J_{t-dt}' * (u0 o iT_{t-dt})
        m = transfo('map', m0, iT, vs);             % <- u0 o iT_{t-dt}
        m = pointwise_multrans_jacobian(iJ, m);     % <- J_{t-dt}' * (u0 o iT_{t-dt})
        m = jacobian_determinant(iJ) .* m;
        
        % Update velocity (eq.16)
        % v_{t+dt} = K u_{t+dt}
        v = diffeo('mom2vel', m, K);
    end
    fprintf('| end\n');
end

%% Helpers: PULLEXP/PUSHEXP

function J = pull_jacobian(J, x, J_vs)
% FORMAT J = pull_jacobian(J, x, J_vs)
% J     - Jacobian tensor volume ([lat_dim vec_dim vec_dim])
% x     - Small deformation ([lat_dim vec_dim])
% J_vs  - Voxel size of the Jacobian tensor lattice
%
% Evaluate the Jacobian tensor at a new position (J o x)
    dim_J   = size(J);
    vec_dim = dim_J(end);
    lat_dim = dim_J(1:vec_dim);
    % Distribute tensor coefficeints along the last dimension
    J = reshape(J, [lat_dim vec_dim * vec_dim]);
    % Interpolate J (1st order, circulant) and evaluate at x
    J = transfo('map', J, x, J_vs, 1, 1);
    J = reshape(J, dim_J);
end

function K = pointwise_multrans_jacobian(J, x)
% FORMAT K = pointwise_multrans_jacobian(J, x)
% J - Jacobian tensor volume ([lat_dim vec_dim vec_dim])
% x - Vector field on the same lattice as J ([lat_dim vec_dim []])
%
% Performs J' * x at each point of the lattice
    dim_J = size(J);
    vec_dim = dim_J(end);
    lat_dim = dim_J(1:vec_dim);
    dim_x = size(x);
    dim_x_point = dim_x((length(lat_dim)+1):end);
    x = reshape(x, [prod(lat_dim) dim_x_point]);
    J = reshape(J, [], vec_dim, vec_dim);
    K = zeros(size(x));
    for i=1:size(x,1)
        jac = reshape(J(i,:,:), [vec_dim vec_dim]);
        p = reshape(x(i,:,:), vec_dim, size(x,3));
        K(i,:,:) = jac' * p;
    end
    K = reshape(K, [lat_dim dim_x_point]);
end

function D = jacobian_determinant(J)
% FORMAT D = jacobian_determinant(J)
% J - Jacobian tensor volume ([lat_dim vec_dim vec_dim])
%
% Returns the jacobian determinant at each point of the lattice.
    dim_J = size(J);
    vec_dim = dim_J(end);
    lat_dim = dim_J(1:vec_dim);
    J = reshape(J, [prod(lat_dim) vec_dim vec_dim]);
    D = zeros([size(J,1) 1]);
    for i=1:size(J,1)
        D(i) = det(reshape(J(i,:,:), [vec_dim vec_dim]));
    end
    D = reshape(D, lat_dim);
end

function J = pointwise_expm(J)
% FORMAT E = pointwise_expm(J)
% J - Tensor field ([lat_dim vec_dim vec_dim])
% E - Pointwise matrix exponential of the tensor field
%
% Pointwise matrix exponential.
    dim_J = size(J);
    vec_dim = dim_J(end);
    lat_dim = dim_J(1:vec_dim);
    J = reshape(J, [prod(lat_dim) vec_dim vec_dim]);
    for i=1:size(J,1)
        J(i,:,:) = expm(reshape(J(i,:,:), [vec_dim vec_dim]));
    end
    J = reshape(J, dim_J);
end

function J = pointwise_invert(J)
% FORMAT iJ = pointwise_expm(J)
% J  - Tensor field ([lat_dim vec_dim vec_dim])
% iJ - Pointwise matrix inversion of the tensor field
%
% Pointwise matrix inversion.
    dim_J = size(J);
    vec_dim = dim_J(end);
    lat_dim = dim_J(1:vec_dim);
    J = reshape(J, [prod(lat_dim) vec_dim vec_dim]);
    for i=1:size(J,1)
        J(i,:,:) = inv(reshape(J(i,:,:), [vec_dim vec_dim]));
    end
    J = reshape(J, dim_J);
end

function out = push_field(T, v, vs)
% FORMAT v = pushfield(T, v)
%
% Push values of v to their (nearest) position on the transformed lattice.
% This is a rough approximation of v(T^{-1}).
    dim = size(T);
    im_dim = dim(end);
    lat_dim = dim(1:im_dim);
    
    % Convert T to voxels (initially in mm)
    for i=1:im_dim
        T = arrayutils('assign_slice', T, ...
                       arrayutils('select_slice', T, i) / vs(i), ...
                       i);
    end
    % Approximate to nearest value
    T = round(T);
    
    % Push values
    dimv = size(v);
    if length(dimv) > im_dim
        vec_dim = dimv(end);
    else
        vec_dim = 1;
    end
    
    v = reshape(v, [prod(lat_dim) vec_dim]);
    T = reshape(T, [prod(lat_dim) im_dim]);
    out = zeros(size(v));
    n = zeros([prod(lat_dim) 1], 'uint16');
    % ^ to keep track of the number of values assigned to a same point of 
    % the output lattice
    
    % Ugly loop over the lattice
    for i=1:prod(lat_dim)
        % convert T(i) to a linear index
        j = arrayutils('select_slice', T, i, 1);
        j = mod(j-1, lat_dim)+1; % circulant boundary
        j = num2cell(j);
        j = sub2ind(lat_dim, j{:});
        out(j,:) = out(j,:) + v(i,:);
        n(j) = n(j) + 1;
    end
    % Normalize
    n = double(n);
    for i=1:vec_dim
        out = arrayutils('assign_slice', out, ...
                         arrayutils('select_slice', out, i) ./ n, ...
                         i);
    end
    clear n
    out(isnan(out)) = 0;
    out = reshape(out, dimv);
end

%% Functions: PENALTY/VEL2MOM

function m = vel2mom(v, L)
% FORMAT m = diffeo('vel2mom', v, L)
% v - Velocity field. It can be provided either in the form of a
%     vector (1D-array) or of a vector field (ND-array).
% L - Regularization matrix. It is the inverse of the prior
%     covariance matrix and can be generated with inv_smooth_prior.
% m - Output momentum field: m = L * v. Same dimensions (1D or ND) as the
%     input.
    m = v;
    dim = size(v);
    if length(dim) > 1
        % Velocity is in volume form
        [m, lat_dim] = arrayutils('field2v', m);
    end
    m = L * m;
    if length(dim) > 1
        % Velocity is in volume form
        m = arrayutils('v2field', m, lat_dim);
    end
end

function L = penalty(lat_dim, lat_vs, param)
% FORMAT L = diffeo('penalty', lat_dim, lat_vs, param)
% lat_dim  - Dimensions of the field lattice ([nx ny ...])
% lat_vs   - Voxel size of the field lattice in mm/vox ([vx vy ...])
% param[1] - Weight of the absolute displacement penalty (mm^{-2})
% param[2] - Lambda parameter of the membrane energy.
% param[3] - Lambda parameter of the bending energy (mm^2).
% param[4] - Mu parameter of the linear-elastic energy.
% param[5] - Lambda parameter of the linear-elastic energy.
%
% Computes the inverse of the covariance matrix of the smooth prior 
% in the LDDMM framwork (i.e., the posdef, self-adjoint, differential 
% operator that defines a metric in the velocity space). It penalizes 
% absolute deformations and membrane, bending and linear-elastic energies. 
% See Ashburner, "A fast diffeomorphic image registration algorithm",
% NeuroImage (2007)
%
% Robustness is improved by averaging sum of square penalties obtained by
% different approximations of the Jacobian (right or left approximation),
% as done in spm_diffeo.
    im_dim = length(lat_dim);
    if nargin < 3
        param = [0.0001 0.001 0.2 0.05 0.2];
    end
    if nargin < 2
        lat_vs = ones(1, im_dim);
    end
    
    lambda_abs      = param(1); % sqrt(lam) ~ 1/mm
    lambda_membrane = param(2); % sqrt(lam) ~ 1
    lambda_bending  = param(3); % sqrt(lam) ~ mm
    mu_elastic      = param(4); % sqrt(mu)  ~ 1
    lambda_elastic  = param(5); % sqrt(lam) ~ 1
    
    % Compute all different combination of Jacobian estimate
    % ([Jx+ Jy+], [Jx+ Jy-], [Jx- Jy+], [Jx- Jy-], ...)
    dirs = combdir(im_dim);

    % Start with a zero matrix
    L = sparse(prod(lat_dim) * im_dim, prod(lat_dim) * im_dim);
    % Absolute displacement in each component (~= euclidian dist)
    if lambda_abs > 0
        E = moperators('euclidian_dist', lat_dim, lat_vs);
        L = L + lambda_abs * (E' * E);
    end
    % Membrane energy
    if lambda_membrane > 0
        for i=1:numel(dirs)
            J = moperators('jacobian', lat_dim, lat_vs, im_dim, dirs{i});
            L = L + lambda_membrane / numel(dirs) * (J' * J);
        end
    end
    % Bending energy
    if lambda_bending > 0
        for i=1:numel(dirs)
            H = moperators('hessian', lat_dim, lat_vs, im_dim, dirs{i});
            L = L + lambda_bending / numel(dirs) * (H' * H);
        end
    end
    % Linear-elastic (symmetric part of the jacobian)
    if mu_elastic > 0
        for i=1:numel(dirs)
            S = moperators('symjac', lat_dim, lat_vs, im_dim, dirs{i});
            L = L + 2 * mu_elastic / numel(dirs) * (S' * S);
        end
    end
    % Linear-elastic (divergence)
    if lambda_elastic > 0
        for i=1:numel(dirs)
            D = moperators('divergence', lat_dim, lat_vs, dirs{i});
            L = L + lambda_elastic / numel(dirs) * (D' * D);
        end
    end
end

%% Helpers: PENALTY/VEL2MOM

function CD = combdir(im_dim)
% Computes all possible combinations of directions for penalty operators
    sets = mat2cell(repmat([-1 1], im_dim, 1), ones(1,im_dim), 2);
    out = cell(im_dim, 1);
    [out{1:im_dim}] = ndgrid(sets{:});
    CD = cell(numel(out{1}), 1);
    for i=1:numel(out{1})
        combination = zeros(1, im_dim);
        for j=1:im_dim
            combination(j) = out{j}(i);
        end
        CD{i} = combination;
    end
end

function F = kernel(L, lat_dim, lat_vs)
% FORMAT K = diffeo('kernel', L, lat_dim)
% FORMAT K = diffeo('kernel', param, lat_dim, lat_vs)
% L         - Regularization matrix (i.e. inverse of the prior covariance)
% param     - It is also possible to provide the penalty parameters instead
%             of the full matrix.
% lat_dim   - Dimension of the field lattice. Necessary because it is not
%             explicitely stored in L and is needed to extract the 
%             convolution operators from L.
% (lat_vs)  - Voxel size of the lattice (default: 1).
% F         - Direct convolution operators, equivalent to L.
    if nargin < 3
        lat_vs = ones(size(lat_dim));
    end
    if issame(size(L), [1 5])
        % Params were provided, not the matrix
        L = penalty(lat_dim, lat_vs, param);
    end
    im_dim = length(lat_dim);
    N = prod(lat_dim);
    F = cell(im_dim, im_dim);
    for i=1:im_dim
        ind_i = (N*(i-1) + 1):(N*i); 
        for j=1:im_dim
            ind_j = (N*(j-1) + 1);
            F{i, j} = reshape(full(L(ind_i, ind_j)), lat_dim); 
        end
    end
end

%% Functions GREENS/MOM2VEL

function K = greens(L, lat_dim, lat_vs)
% FORMAT K = diffeo('greens', L, lat_dim)
% FORMAT K = diffeo('greens', param, lat_dim, lat_vs)
% L         - Regularization matrix (i.e. inverse of the prior covariance)
% param     - It is also possible to provide the penalty parameters instead
%             of the full matrix.
% lat_dim   - Dimension of the field lattice. Necessary because it is not
%             explicitely stored in L and is needed to extract the 
%             convolution operators from L.
% (lat_vs)  - Voxel size of the lattice. Only needed when parameters - and
%             not the full matrix - is provided. (default: 1)
% K         - Inverse of L, obtained and stored as convolution operators in 
%             Fourier domain.
%
% Computes the inverse K of L in Fourier space.
% Gets the direct convolution operator for each component and transform
% it to the Fourier domain (since m_i = sum_i conv(c_{i,j}, v_j), 
% where m_i is the i-th momentum component and v_j is the j-th velocity 
% component). Then, makes use of the block structure of the kernel to
% invert it.
    if nargin < 3
        K = kernel(L, lat_dim);
    else
        K = kernel(L, lat_dim, lat_vs);
    end
    im_dim = length(lat_dim);
    for i=1:im_dim
        for j=1:im_dim
            K{i, j} = real(fftn(K{i, j})); 
        end
    end
    K = invert_direct_kernel(K);
end

function v = mom2vel(m, K, lat_dim)
% FORMAT v = diffeo('mom2vel', m, K, lat_dim)
% m         - Momentum field. It can be provided either in the form of a
%             vector (1D-array, providing lat_dim is then mandatory) or of
%             a vector field (ND-array).
% K         - L's inverse (or Green's function). It is stored in the form
%             of a set of convolution kernels in Fourier space and can be
%             generated with greens_kernel.
% (lat_dim) - Dimensions of the field lattice.
% v         - Velocity field. Same dimensions (1D or ND) as the input.
    v = m;
    if nargin == 3
        % Moment is in vector form
        im_dim = length(vel_dim);
        v = reshape(v, lat_dim, []);
    end
    % 1. ND Fourier transform of the momentum
    v = fft_vector_field(v);
    % 2. Inverse of the convolutions in Fourier domain
    v = inverse_convolution_fourier(v, K);
    % 3. Inverse ND Fourier transform of the velocity
    v = ifft_vector_field(v);
    if nargin == 3
        % Moment is in vector form
        v = reshape(v, [], im_dim);
    end
end

%% Helpers: GREENS/MOM2VEL

function f = fft_vector_field(v)
    % Performs ND Fourier transform on each component of the vector field
    dim = size(v);
    vec_dim = dim(end);
    f = v;
    for i=1:vec_dim
        f = arrayutils('assign_slice', f, fftn(arrayutils('select_slice', f, i)), i);
    end
end

function v = ifft_vector_field(f)
    % Performs ND inverse Fourier transform on each component of the vector 
    % field (supposded symmetric)
    dim = size(f);
    vec_dim = dim(end);
    v = f;
    for i=1:vec_dim
        v = arrayutils('assign_slice', v, ifftn(arrayutils('select_slice', v, i), 'symmetric'), i);
    end
end

function fv = inverse_convolution_fourier(fm, K)
% FORMAT fv = inverse_convolution_fourier(fm, fK)
% fm - Fourier transform of the momentum (dimensions [lat_dim vec_dim])
% K  - Fourier transform of the inverse convolution kernels, can be
%      generated with greens_kernel.
% fv - Fourier transform of the velocity field.
%
% Compute Fv (velocity in Fourier space) as Fv_i = sum_j K{i,j} Fm_j
    dim = size(fm);
    lat_dim = dim(1:end-1);
    vec_dim = dim(end);
    fv = zeros([lat_dim vec_dim]);
    % Apply kernel
    for i=1:vec_dim
        for j=1:vec_dim
            convolved = arrayutils('select_slice', fm, j) .* K{i,j};
            fv = arrayutils('assign_slice', fv, arrayutils('select_slice', fv, i) + convolved, i);
        end
    end
end

function [K, check] = invert_direct_kernel(F)
% FORMAT K = invert_direct_kernel(iK)
% F     - Cell array containing the direct convolution operators of L in 
%         the Fourier domain.
% check - Gauss-Jordan transform of the input, should be a block-identity 
%         matrix if all went well.
%
% Make use of the blocks-of-diagonals structure of iK to invert it with 
% block-Gauss-Jordan elimination.
    N = size(F, 1);
    K = cell(N);
    % Initialize the block identity matrix
    for i=1:N
        for j=1:N
            if i == j
                K{i,j} = ones(size(F{i,j}));
            else
                K{i,j} = zeros(size(F{i,j}));
            end
        end
    end
    % Gauss-Jordan
    for i=1:N
        factor = F{i,i};
        for j=1:N
            K{i,j} = K{i,j} ./ factor;
            F{i,j} = F{i,j} ./ factor;
        end
        for k=1:N
            if k ~= i
                factor = F{k,i};
                for j=1:N
                    K{k,j} = K{k,j} - K{i,j} .* factor;
                    F{k,j} = F{k,j} - F{i,j} .* factor;
                end
            end
        end
    end
    check = F;
end