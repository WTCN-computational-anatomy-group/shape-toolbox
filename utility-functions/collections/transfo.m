function varargout = transfo(id, varargin)
% FORMAT id  = transfo('idmap', lat_dim, (lat_vs))
% FORMAT t   = transfo('trl2map', T, lat_dim, (lat_vs))
% FORMAT t   = transfo('aff2map', A, lat_dim, (lat_vs))
% FORMAT y   = transfo('map', T, x, (interp))
% FORMAT t   = transfo('compose', t1, (vs1), ..., tn, (vsn), ...)
% FORMAT out = transfo('warp', ima, T, ...)
% FORMAT ok  = transfo('test')
%
% Collection of tools for manipulating transformation maps.
%
% FORMAT help transfo>function
% Returns the help file of the selected function.
    if nargin == 0
        id = 'test';
    end
    switch id
        case 'idmap'
            [varargout{1:nargout}] = idmap(varargin{:});
        case 'trl2map'
            [varargout{1:nargout}] = trl2map(varargin{:});
        case 'aff2map'
            [varargout{1:nargout}] = aff2map(varargin{:});
        case 'map'
            [varargout{1:nargout}] = map(varargin{:});
        case 'compose'
            [varargout{1:nargout}] = compose(varargin{:});
        case 'warp'
            [varargout{1:nargout}] = warp(varargin{:});
        case 'test'
            [varargout{1:nargout}] = test(varargin{:});
        case 'help'
            [varargout{1:nargout}] = help_func(varargin{:});
    end
end

%% Functions

function id = idmap(lat_dim, lat_vs)
% FORMAT id = transfo('idmap', lat_dim, lat_vs)
% lat_dim - Dimensions of the lattice on which to compute the map
% lat_vs  - Voxel size of the lattice
%
% Returns the identity transform on a lattice.
    if nargin < 2
        lat_vs = ones(size(lat_dim));
    end
    im_dim = length(lat_dim);
    id = zeros([lat_dim im_dim]);
    for i=1:im_dim
        tmp_dim = lat_dim;
        tmp_dim(i) = 1;
        one_dim = ones(1, length(lat_dim));
        one_dim(i) = lat_dim(i);
        id = arrayutils('assign_slice', id, ...
                        repmat(reshape((1:lat_dim(i)) * lat_vs(i), one_dim), tmp_dim), ...
                        i);
    end
end

function y = map(T, x, T_vs, interp, circulant)
% FORMAT y = transfo('map', T, x, interp)
% T         - Discretized transformation on a lattice of dimensions 
%             [nx ny ...]
% x         - Array of coordinates of dimensions [M N ... D] with D the 
%             dimension of the transformation (2D/3D...)
% (T_vs)    - Voxel size of the transformation lattice (default: 1).
% (interp ) - Interpolation degree of the transformation map (default: 1).
% y         - Transformed coordinates
%
% Apply a tranformation map to a coordinate, a set of coordinates or
% another transformation map (i.e., compose the maps).
    if nargin < 3
        T_vs = 1;
    end
    if nargin < 4
        interp = 1;
    end
    if nargin < 5
        circulant = 1;
    end
    % Reshape x
    dim_x = size(x);
    ndim_x = length(dim_x);
    dimvec_x = dim_x(ndim_x);
    dimlat_x = dim_x(1:dimvec_x);
    dimlat_x_3d = dimlat_x;
    dim_T = size(T);
    dimvec_T = dim_T(end);
    dimlat_T = dim_T(1:dimvec_x);
    switch ndim_x
        case 1
            x = reshape(x, [1 1 1 dimvec_x]);
            dimlat_x_3d = [1 1 1];
        case 2
            x = reshape(x, [size(x, 1) 1 1 dimvec_x]);
            dimlat_x_3d = [size(x, 1) 1 1];
        case 3
            x = reshape(x, [size(x, 1) size(x, 2) 1 dimvec_x]);
            dimlat_x_3d = [size(x, 1) size(x, 2) 1];
    end
    if numel(T_vs) < dimvec_x
        T_vs = padarray(T_vs, [0 dimvec_x - numel(T_vs)], 1, 'post');
    end
    % Convert coordinates from mm to voxels
    for i=1:dimvec_x
        x(:,:,:,i) = x(:,:,:,i) / T_vs(i);
    end
    % Generate spline coefficients
    if numel(interp) < 3
        interp = padarray(interp, [0 3 - numel(interp)], 'replicate', 'post');
    end
    circulant = double(circulant);
    if numel(circulant) < 3
        circulant = padarray(circulant, [0 3 - numel(circulant)], 'replicate', 'post');
    end
    y = zeros([dimlat_x_3d dimvec_T]);
    % Interpolate the transformation (T) at the provided coordinates (x)
    for i=1:dimvec_T
        if issame(interp, zeros([1 3]))
            % If interpolation order == [0 0 0], I perform nearest neighbour myself
            % to avoid bad precision issueswith bsplins
            nrst_x = floor(x + 0.5);
            nrst_y = arrayutils('select_slice', T, i);
            nrst_x = reshape(nrst_x, [], dimvec_x);
            for j=1:dimvec_x
                nrst_x(:,j) = mod(nrst_x(:,j)-1, dimlat_T(j))+1;
            end
            switch dimvec_x
                case 2
                    nrst_x = sub2ind(size(nrst_y), nrst_x(:,1), nrst_x(:,2));
                case 3
                    nrst_x = sub2ind(size(nrst_y), nrst_x(:,1), nrst_x(:,2), nrst_x(:,3));
            end
            nrst_y = nrst_y(nrst_x);
            nrst_y = reshape(nrst_y, dim_x(1:(end-1)));
            y = arrayutils('assign_slice', y, nrst_y, i);
        else
            coeffs = spm_diffeo('bsplinc', ...
                                single(arrayutils('select_slice', T, i)), ...
                                [interp circulant]);
            one_y = spm_diffeo('bsplins', coeffs, single(x), [interp circulant]);
            one_y(isnan(one_y)) = 0;
            y = arrayutils('assign_slice', y, one_y, i);
        end
    end
    % Reshape y
    y = reshape(y, [dimlat_x dimvec_T]);
end

function t = compose(varargin)
% FORMAT t = transfo('compose', t1, (vs1), t2, (vs2), ..., tn, (vsn),
%                               (interp))
% t1, ..., tn     - n transformation maps.
% (vs1, ..., vsn) - voxel sizes of the transformations.
% (interp)        - Degree to use to interpolate maps (default: 1).
% t               - Composed transformation (t1 o t2 o ... o tn).
%
% Compose transformation maps from right to left. The output voxel size
% will be that of the left most transformation.
    % If the last argument is a scalar, it is the interpolation degree
    t_dim = size(varargin{1,1});
    t_dim = t_dim(end);
    vs_provided = numel(varargin{1,2}) <= t_dim;
    ntransfo = nargin;
    if numel(varargin{1,end}) <= t_dim 
        if ~vs_provided || numel(varargin{1,end-1}) <= t_dim
            d = varargin{1,end};
            varargin = {varargin{1,1:(end-1)}};
            ntransfo = ntransfo - 1;
        else
            d = 1;
        end
    else
        d = 1;
    end
    if vs_provided
        ntransfo = ntransfo/2;
    end
    function j = gettransfo(i)
        j = ntransfo - i + 1;
        if vs_provided
            j = 2*j - 1;
        end
    end
    % Compose right to left
    t = varargin{1,gettransfo(1)};
    for i=2:ntransfo
        T = varargin{1,gettransfo(i)};
        dim_T = size(T);
        dimvec_T = dim_T(end);
        dimlat_T = dim_T(1:dimvec_T);
        if vs_provided
            id = idmap(dimlat_T, varargin{1,gettransfo(i)+1});
            T = T - id;
            t = map(T, t, varargin{1,gettransfo(i)+1}, d) + t;
        else
            id = idmap(dimlat_T);
            T = T - id;
            t = map(T, t, d) + t;
        end
            
    end
end

function out = warp(ima, T, varargin)
% FORMAT out = transfo('warp', ima, T,
%                              'out_dim',    out_dim,
%                              'out_vs',     out_vs,
%                              'interp_ima', interp_ima,
%                              'interp_T',   interp_T,
%                              'T_vs',       T_vs,
%                              'ima_vs',     ima_vs)
% ima           - ND-array of dimensions [nx ny ...] (i.e., a 
%                 scalar image) or [nx ny ... k] (i.e., a vector image).
% T             - Transformation map of dimensions ([nx ny ... d]). It must
%                 map output coordinates towards image coordinates (i.e., 
%                 the invert of what's natural).
% (out_dim)     - Dimensions of the output lattice (default: same as T).
% (out_vs)      - Voxel size of the output lattice (default: same as T).
% (interp_ima)  - Order of the image interpolating spline (default: 1).
% (interp_T)    - Order of the transformation interpolating spline
%                 (default: 1).
% (T_vs)        - Voxel size of the interpolating spline.
% (ima_vs)      - Voxel size of the input image.
%
% Warp an image towards an output lattice.
    T_dim = size(T);
    ima_dim = size(ima);
    vol_dim = T_dim(end);
    % Default arguments
    opt             = struct;
    opt.out_dim     = [];
    opt.out_vs      = [];
    opt.interp_ima  = 1;
    opt.interp_T    = 1;
    opt.T_vs        = 1;
    opt.ima_vs      = 1;
    % Parse arguments
    opt             = parse_varargin(varargin, opt);
    if isempty(opt.out_dim)
        opt.out_dim = T_dim(1:(end-1));
    end
    % Select scalar or vector case
    is_vector_image = length(ima_dim) == vol_dim + 1;
    if is_vector_image % vector case
        vec_dim = ima_dim(end);
        out = zeros([opt.out_dim vec_dim]);
        for i=1:vec_dim
            slice = warp_scalar_image(arrayutils('select_slice', ima, i), T, varargin{:});
            out = arrayutils('assign_slice', out, slice, i);
        end
    else % scalar case
        out = warp_scalar_image(ima, T, varargin{:});
    end 
end

function out = warp_scalar_image(ima, T, varargin)
% FORMAT out = transform_scalar_image(ima, T,
%                                     'out_dim',    out_dim,
%                                     'out_vs',     out_vs,
%                                     'interp_ima', interp_ima,
%                                     'interp_T',   interp_T,
%                                     'T_vs',       T_vs,
%                                     'ima_vs',     ima_vs)
% ima           - Scalar ND-array of dimensions [nx ny ...] (i.e., a 
%                 scalar image).
% T             - Transformation map of dimensions ([nx ny ... d]). It must
%                 map output coordinates towards image coordinates (i.e., 
%                 the invert of what's natural).
% (out_dim)     - Dimensions of the output lattice (default: same as T).
% (out_vs)      - Voxel size of the output lattice (default: same as T).
% (interp_ima)  - Order of the image interpolating spline (default: 1).
% (interp_T)    - Order of the transformation interpolating spline
%                 (default: 1).
% (T_vs)        - Voxel size of the interpolating spline.
% (ima_vs)      - Voxel size of the input image.
%
% Warp a scalar image towards an output lattice.
    % Default arguments
    opt             = struct;
    opt.out_dim     = [];
    opt.out_vs      = [];
    opt.interp_ima  = 1;
    opt.interp_T    = 1;
    opt.T_vs        = 1;
    opt.ima_vs      = 1;
    % Parse arguments
    opt             = parse_varargin(varargin, opt);
    % Process arguments
    T_dim = size(T);
    if numel(opt.interp_ima) < 3
        opt.interp_ima = padarray(opt.interp_ima, ...
                                  [0 3 - numel(opt.interp_ima)], ...
                                  'replicate', 'post');
    end
    if numel(opt.interp_T) < 3
        opt.interp_t = padarray(opt.interp_T, ...
                                [0 3 - numel(opt.interp_T)], ...
                                'replicate', 'post');
    end
    boundaries = ones(1, 3);
    if numel(opt.T_vs) < T_dim(end)
        opt.T_vs = padarray(opt.T_vs, ...
                            [0 T_dim(end) - numel(opt.T_vs)], ...
                            'replicate', 'post');
    end
    if numel(opt.ima_vs) < T_dim(end)
        opt.ima_vs = padarray(optima_vsT_vs, ...
                              [0 T_dim(end) - numel(opt.ima_vs)], ...
                              'replicate', 'post');
    end
    if isempty(opt.out_dim)
        opt.out_dim = T_dim(1:(end-1));
    end
    if isempty(opt.out_vs)
        opt.out_vs = opt.T_vs;
    end
    
    % Image spline transform
    ima_coeff = spm_diffeo('bsplinc', single(ima), ...
                           [opt.interp_ima boundaries]);
    % Coordinates
    out_coordinates = idmap(opt.out_dim, opt.out_vs);
    out_coordinates = compose(out_coordinates, opt.out_vs, ...
                                              T, opt.T_vs, opt.interp_T);
    % Reshape coordinates
    switch T_dim(end)
        case 1
            out_coordinates = reshape(out_coordinates, [size(out_coordinates, 1) 1 1 T_dim(end)]);
        case 2
            out_coordinates = reshape(out_coordinates, [size(out_coordinates, 1) size(out_coordinates, 2) 1 T_dim(end)]);
    end
    % Convert coordinates from mm to voxels
    for i=1:T_dim(end)
        out_coordinates(:,:,:,i) = out_coordinates(:,:,:,i) / opt.ima_vs(i);
    end
    % Interpolate image on output grid
    out = spm_diffeo('bsplins', ima_coeff, single(out_coordinates), ...
                     [opt.interp_ima boundaries]);
end

function t = trl2map(T, lat_dim, lat_vs)
% FORMAT t = transfo('trl2map', T, lat_dim, lat_vs)
% T         - Array of translation parameters
% lat_dim   - Dimensions of the map lattice
% lat_vs    - Voxel size of the map lattice
%
% Generate a transformation map from a translation.
    D = length(T);
    A = [eye(D) T(:)];
    t = aff2map(A, lat_dim, lat_vs);
end

function t = aff2map(A, lat_dim, lat_vs)
% FORMAT t = transfo('aff2map', A, lat_dim, lat_vs)
% A         - Affine transformation matrix ([D * D+1] or [D+1 * D+1])
% lat_dim   - Dimensions of the map lattice
% lat_vs    - Voxel size of the map lattice
%
% Generate a transformation map from an affine transformation.
    t = idmap(lat_dim, lat_vs);
    t = [reshape(t, prod(lat_dim), length(lat_dim)) ones(prod(lat_dim), 1)];
    [M,N] = size(A);
    if M == N-1
        A = [A ; zeros(1, length(lat_dim) + 1)];
    end
    [M,N] = size(A);
    if M ~= length(lat_dim) + 1 || N ~= length(lat_dim) + 1
        error('A should be a %iD affine matrix in homogenate coordinates', length(lat_dim));
    end
    t = (A * t')';
    t = t(:,1:(end-1));
    t = reshape(t, [lat_dim length(lat_dim)]);
end

%% Test
function ok = test
% FORMAT ok = transfo('test')
%
% Launches unit tests. Should return true.
ok = true;

dim = [28 28];
vs = [1 1];

id = idmap(dim, vs);
tr2 = id + padarray(2 * ones(dim), [0 0 1], 0, 'post'); % 2mm x-translation
tr3 = id + padarray(3 * ones(dim), [0 0 1], 0, 'post'); % 3mm x-translation
x = [3 5; 7 10];

y0 = map(id, x, vs);
ok = ok && isequal(x, y0);
y2 = map(tr2, x, vs);
y3 = map(tr3, x, vs);
tr5 = compose(tr2, vs, tr3, vs); % should be 2+3=5mm x-translation
y5 = map(tr5, x, vs);
yt = map(tr2, y3, vs);
ok = ok && isequal(y5, yt);

% Generate random image
vel_dim = [28 28];
vel_vs = [1 1];
im_dim = 2;
v = zeros([vel_dim, im_dim]);
vel_fwhm = 20;
for i=1:im_dim
    if im_dim == 2
        v(:,:,i) = spm_conv(randn(vel_dim), vel_fwhm);
    elseif im_dim == 3
        tmp = zeros(vel_dim);
        spm_smooth(randn(vel_dim), tmp, vel_fwhm);
        v(:,:,:,i) = tmp;
    end
end

% Warp image
A = [ .3 0 -2 ;
      0 1 0 ];
trn = aff2map(A, vel_dim, vel_vs);
v2 = warp(v, trn, 'ima_vs', vel_vs, 'T_vs', vs, 'out_vs', vel_vs, 'out_dim', vel_dim);
v3 = warp(v, id, 'ima_vs', vel_vs, 'T_vs', vs, 'out_vs', vel_vs, 'out_dim', vel_dim);
end