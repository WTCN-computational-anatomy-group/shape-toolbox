function varargout = moperators(id, varargin)
% FORMAT J = moperators('jacobian', lat_dim, lat_vs, vec_dim)
% FORMAT H = moperators('hessian', lat_dim, lat_vs, vec_dim)
% FORMAT H = moperators('laplacian', lat_dim, lat_vs)
% FORMAT D = moperators('euclidian_dist', lat_dim, lat_vs)
% FORMAT A = moperators('divergence', lat_dim, lat_vs)
% FORMAT S = moperators('symjac', lat_dim, lat_vs, vec_dim)
% 
% Collection of matricial operators. For exemple, let v be the vectorized 
% version of a scalar or vector field, and J be the jacobian operator, then
% J * v returns a vectorized version of the Jacobian matrix of v.
% 
% FORMAT help moperators>function
% Returns the help file of the selected function.
    switch id
        case 'jacobian'
            [varargout{1:nargout}] = jacobian(varargin{:});
        case 'hessian'
            [varargout{1:nargout}] = hessian(varargin{:});
        case 'euclidian_dist'
            [varargout{1:nargout}] = euclidian_dist(varargin{:});
        case 'laplacian'
            [varargout{1:nargout}] = laplacian(varargin{:});
        case 'divergence'
            [varargout{1:nargout}] = divergence(varargin{:});
        case 'symjac'
            [varargout{1:nargout}] = symjac(varargin{:});
        otherwise
            error('No function named %s\n', string(id));
    end
end


%% Matricial operators

function J = jacobian(lat_dim, lat_vs, vec_dim, dir)
% FORMAT J = moperators('jacobian', lat_dim, lat_vs, vec_dim, dir)
% lat_dim   - Dimensions of the lattice ([nx ny ...]).
% (lat_vs)  - Voxel size of the lattice (default: 1).
% (vec_dim) - Dimension of the vector space (default: 1).
% (dir)     - Approximation used in each direction ([dx dy ...]):
%             * -1: Left-hand approximation
%             *  1: Right-hand approximation
%             *  0: Mean approximation (default)
% J         - Jacobian matrix operator s.t. Jac(v) = J * v
%
% Return a matrix that allows estimating the Jacobian of a vector
% field at a given point. The partial derivative at a point is 
% estimated by:
% > dv/dx (p) ~= ( v(p + dx) - v(p) ) / dx , 
% where dx is one voxel.
% The returned matrix J can be used through 
% > Jac(v) = reshape(J * v(:), [lat_dim field_dim vec_dim]),
% where Jac(v)(x) is the Jacobian at point x.
% If vector_dim is provided, the field v is considered to be vectorized
% over both the lattice and the vectors (v = [v1; v2; v3]), and the 
% returned matrix J has dimensions 
% (prod(lattice_dim) * vector_dim) x (prod(lattice_dim)).
% Else, it is considered to be vectorized only over the lattice 
% (v = [v1 v2 v3]), and the returned matrix J has dimensions 
% (prod(lattice_dim)) x (prod(lattice_dim)).
%
% The sum of square derivatives at point v is: v'*(J'*J)*v
    
    if nargin < 4
        dir = zeros(1, length(lat_dim));
    end
    if nargin < 3
        vec_dim = 1;
    end
    if nargin < 2
        lat_vs = ones(1, length(lat_dim));
    end
        
    % Multidimensional gradient functor
    CG = multi_gradient(lat_dim, dir, lat_vs);
    
    % Concatenate gradient functors
    % > CG * f allows to compute the gradient at f, where f is a scalar 
    %   field.
    CG = cell2mat(CG); 
    
    % Vector gradient functor
    % > J * v allows to compute the jacobian at v, where v is a vector 
    %   field.
    J = kron(speye(vec_dim), CG);
end

function H = hessian(lat_dim, lat_vs, vec_dim, dir)
% FORMAT H = moperators('hessian', lat_dim, lat_vs, vec_dim)
% lat_dim   - Dimensions of the lattice ([nx ny ...]).
% (lat_vs)  - Voxel size of the lattice (default: 1).
% (vec_dim) - Dimension of the vector space (default: 1).
% (dir)     - Approximation used in each direction ([dx dy ...]):
%             * -1: Left-hand approximation
%             *  1: Right-hand approximation
%             *  0: Mean approximation (default)
% H         - Hessian matrix operator s.t. Hess(v) = H * v
%
% Returns a matrix that allows estimating the Hessian of a vector 
% field at a given point. The partial derivative at a point is 
% estimated by d2v/dxdy (p) ~= ( dv/dx(p + dy) - dv/dx(p) ) / dy , 
% where dy is one voxel.
% The returned matrix J can be used through 
% > Hess(v) = reshape(H * v(:), [lat_dim n_comb vec_dim]), 
% where Hess(v)(i,:,:,j) is the Hessian at point v_i (stored sparse) 
% for component j, and n_comb = lat_dim * (lat_dim+1) / 2.
% If vector_dim is provided, the field v is considered to be vectorized
% over both the lattice and the vectors (v = [v1; v2; v3]), and the 
% returned matrix H has dimensions 
% (prod(lattice_dim) * vector_dim) x (prod(lattice_dim)).
% Else, it is considered to be vectorized only over the lattice 
% (v = [v1 v2 v3]), and the returned matrix J has dimensions 
% (prod(lattice_dim)) x (prod(lattice_dim)).
%
% The sum of square second derivatives at point v is: v'*(H'*H)*v
    
    if nargin < 4
        dir = zeros(1, length(lat_dim));
    end
    if nargin < 3
        vec_dim = 1;
    end
    if nargin < 2
        lat_vs = ones(1, length(lat_dim));
    end
    im_dim = length(lat_dim);
    
    % Multidimensional gradient functor
    CG = multi_gradient(lat_dim, dir, lat_vs);
    
    % Multidimensional hessian functor
    % > CGG{index(i,j)} * f allows to compute the second derivative along 
    %   directions i and j, where f is a scalar field.
    % I don't use the symmetric approximation because the sqrt(2)
    % multiplication causes imprecision errors.
    CGG = cell(im_dim * im_dim, 1);
    index = 1;
    for i=1:im_dim
        for j=1:im_dim
            CGG{index} = CG{i} * CG{j};
            index = index + 1;
        end
    end
    clear CG
    
    % Concatenate hessian functors
    % > CG * f allows to compute the hessian at f, where 
    %   f is a scalar field.
    CGG = cell2mat(CGG); 
    
    % Vector hessian functor
    % > H * v allows to compute the hessian at v, where v is a vector 
    %   field.
    H = kron(speye(vec_dim), CGG);
end

function L = laplacian(lat_dim, lat_vs, dir)
% FORMAT L = moperators('laplacian', lat_dim, lat_vs)
% lat_dim   - Dimensions of the lattice ([nx ny ...]).
% (lat_vs)  - Voxel size of the lattice (default: 1).
% (dir)     - Approximation used in each direction ([dx dy ...]):
%             * -1: Left-hand approximation
%             *  1: Right-hand approximation
%             *  0: Mean approximation (default)
% L         - Laplacian matrix operator s.t. Lap(v) = L * f
%
% Returns a matrix that allows estimating the laplacian of a scalar field
% at a given point. The partial derivative at a point is 
% estimated by d2f/d2x (p) ~= ( df/dx(p + dx) - df/dx(p) ) / dx , 
% where dx is one voxel.
% The returned matrix can be used through
% > Lapl(f) = reshape(L * f(:), [lat_dim]),
% where Lapl(f)(x) is the value of the Laplacian at point x.
%
% The sum of square laplacians at point f is: f'*(L'*L)*f

    if nargin < 3
        dir = zeros(1, length(lat_dim));
    end
    if nargin < 2
        lat_vs = ones(1, length(lat_dim));
    end
    im_dim = length(lat_dim);

    % Multidimensional gradient functor
    CG = multi_gradient(lat_dim, dir, lat_vs);
    
    % Multidimensional diagonal of hessian functor
    % > CGG{i} * f allows to compute the unmixed second derivative along 
    %   direction i,  where f is a scalar field.
    CGG = cell(1, im_dim);
    for i=1:im_dim
        CGG{i} = CG{i} * CG{i};
    end
    clear CG
    
    % Concatenate diagonal of hessian functors
    % > L * f allows to compute the laplacian at f, where 
    %   f is a scalar field.
    L = cell2mat(CGG); 
end

function D = euclidian_dist(lat_dim, lat_vs)
% FORMAT D = moperators('euclidian_dist', lat_dim, lat_vs)
% lat_dim   - Dimensions of the lattice ([nx ny ...]).
% (lat_vs)  - Voxel size of the lattice (default: 1).
% D         - Distance matrix operator s.t. Dist(v) = D * v
%
% Returns a matrix that allows *estimating* the euclidian norm of a  
% vector field at a given point. 
% The returned matrix D can be used through 
% > Dist(v) = reshape(D * v(:), [lat_dim vec_dim]), 
% where Dist(v)([x,j]) is the value (in mm) of component j at point x.
%
% The sum of square distance components at point v is: v'*(D'*D)*v,
% which is not the sum of square euclidian distances!
    
    if nargin < 2
        lat_vs = ones(1, length(lat_dim));
    end
    vector_dim = length(lat_dim);
    
    D = cell(1, vector_dim);
    for i=1:vector_dim
        D{i} = speye(prod(lat_dim));
        if lat_vs(i) ~= 1
            D{i} = D{i} * lat_vs(i);
        end
    end
    D = blkdiag(D{:});
end

function A = divergence(lat_dim, lat_vs, dir)
% FORMAT A = moperators('divergence', lat_dim, lat_vs, dir)
% lat_dim   - Dimensions of the lattice ([nx ny ...]).
% (lat_vs)  - Voxel size of the lattice (default: 1).
% (dir)     - Approximation used in each direction ([dx dy ...]):
%             * -1: Left-hand approximation
%             *  1: Right-hand approximation
%             *  0: Mean approximation (default)
% A         - Divergence matrix operator s.t. Div(v) = A * v
%
% Returns a matrix that allows estimating the divergence of a vector 
% field at a given point. 
% The returned matrix A can be used through 
% > Div(v) = reshape(A * v(:), lat_dim), 
% where Div(v)(x) is the divergence at point x.
%
% The sum of square divergences at point v is: v'*(A'*A)*v

    if nargin < 3
        dir = zeros(1, length(lat_dim));
    end
    if nargin < 2
        lat_vs = ones(1, length(lat_dim));
    end
    
    % Multidimensional gradient functor (column shape)
    CG = multi_gradient(lat_dim, dir, lat_vs);
    
    % Divergence functor (row shape)
    % > A * v allows to compute the divergence at v, where v is a vector 
    %   field.
    A = cell2mat(CG');
end

function S = symjac(lat_dim, lat_vs, vec_dim, dir)
% FORMAT S = moperators('symjac', lat_dim, lat_vs, vec_dim)
% lat_dim   - Dimensions of the lattice ([nx ny ...]).
% (lat_vs)  - Voxel size of the lattice (default: 1).
% (vec_dim) - Dimension of the vector space (default: 1).
% S         - Jacobian matrix operator s.t. SymJac(v) = S * v
%
% Returns a matrix that allows estimating the symmetric part of the
% jacobian of a vector field at a given point. 
% The returned matrix S can be used through 
% > SymJac(v) = reshape(S * v, [lat_dim field_dim vec_dim]),
% where SymJac(v)(i) is the symmetric part of the Jacobian at point 
% v_i.
%
% The sum of square symmetric parts at point v is: v'*(S'*S)*v
    
    im_dim = length(lat_dim);
    J = jacobian(lat_dim, lat_vs, im_dim, dir);
    
    % Compute a matrix that links a the element of a Jacobian (x_i, f_j) 
    % with itself and its symmetric component (x_j, f_i).
    link_sym = spalloc(im_dim * vec_dim, im_dim * vec_dim, im_dim * vec_dim);
    for i=1:im_dim
        for j=1:vec_dim
            if i == j
                link_sym(sub2ind([im_dim vec_dim], i, j), sub2ind([im_dim vec_dim], i, j)) = 1;
            else
                link_sym(sub2ind([im_dim vec_dim], i, j), sub2ind([im_dim vec_dim], i, j)) = 0.5;
                link_sym(sub2ind([im_dim vec_dim], i, j), sub2ind([im_dim vec_dim], j, i)) = 0.5;
            end
        end
    end
    
    % Expand so that it works for any element of the vector field
    S = kron(link_sym, speye(prod(lat_dim)));
    
    % Compute the symmetric part
    S = S * J;
end

%% HELPERS

function G = simple_gradient(dim, dir, vs)
% Simple gradient functor
% > G * f allows to compute the gradient along the first direction at 
%   f, where f is a scalar field.
    if nargin < 3
        vs = 1;
    end
    if nargin < 2
        dir = 0;
    end
    
    if dir > 0
        % Right-hand approximation
        G = spdiags(repmat([-1, 1], dim, 1), [0, 1], dim, dim);
        G(end, 1) = 1; % Ciculating boundaries
    elseif dir < 0
        % Left-hand approximation
        G = spdiags(repmat([-1, 1], dim, 1), [-1, 0], dim, dim);
        G(1, end) = -1; % Ciculating boundaries
    else
        % Mean
        Gp = spdiags(repmat([-1, 1], dim, 1), [0, 1], dim, dim);
        Gp(end, 1) = 1; % Ciculating boundaries
        Gm = spdiags(repmat([-1, 1], dim, 1), [-1, 0], dim, dim);
        Gm(1, end) = -1; % Ciculating boundaries
        G = 0.5 * (Gp + Gm);
    end
    if vs ~= 1
        G = G / vs;
    end
end

function CG = multi_gradient(dim, dir, vs)
% Multidimensional gradient functor
    % > CG{i} * f allows to compute the gradient along direction i at f,  
    %   where f is a scalar field.
    if nargin < 3
        vs = ones(1, length(lat_dim));
    end
    if nargin < 2
        dir = zeros(1, length(lat_dim));
    end
    im_dim = length(dim);
    
    % Create gradient operators for each dimension
    G = cell(im_dim, 1);
    for i=1:im_dim
        G{i} = simple_gradient(dim(i), dir(i), vs(i));
    end
    
    % Extend matrix
    CG = cell(im_dim, 1); 
    for i=1:im_dim
        GorI = cell(im_dim, 1);
        for j=1:im_dim
            if i ~= j
                GorI{j} = speye(dim(j));
            else
                GorI{j} = G{i};
            end
        end
        CG{i} = speye(1);
        for j = 1:im_dim
            CG{i} = kron(GorI{j}, CG{i});
        end
    end
end