function varargout = arrayutils(id, varargin)
% FORMAT vol        = arrayutils('v2vol', v, lat_dim)
% FORMAT [vol, ...] = arrayutils('vol2v', v)
% FORMAT vol        = arrayutils('v2field', v, lat_dim)
% FORMAT [vol, ...] = arrayutils('field2v', v)
% FORMAT out        = arrayutils('select_slice', in, indices, (dim))
% FORMAT out        = arrayutils('assign_slice', in, slice, indices, (dim))
% FORMAT B          = arrayutils('repmat_along', A, dim, n)
%
% Collection of tools for manipulating arrays, volumes, scalar or vector
% fields, etc.
%
% FORMAT help arrayutils>function
% Returns the help file of the selected function.
    switch id
        case 'v2vol'
            [varargout{1:nargout}] = v2vol(varargin{:});
        case 'vol2v'
            [varargout{1:nargout}] = vol2v(varargin{:});
        case 'v2field'
            [varargout{1:nargout}] = v2field(varargin{:});
        case 'field2v'
            [varargout{1:nargout}] = field2v(varargin{:});
        case 'select_slice'
            [varargout{1:nargout}] = select_slice(varargin{:});
        case 'assign_slice'
            [varargout{1:nargout}] = assign_slice(varargin{:});
        case 'repmat_along'
            [varargout{1:nargout}] = repmat_along(varargin{:});
    end
end

function vol = v2vol(v, lat_dim)
% FORMAT vol = arrayutils('v2vol', v, lat_dim)
% v         - 1D array (i.e., a vector)
% lat_dim   - Dimensions of the lattice ([nx ny ...])
% vol       - ND-array of dimensions lat_dim (image or scalar field)
%
% Reshape a vector into a ND-volume (image or scalar field).
    vol = reshape(v, lat_dim);
end

function [vol, lat_dim] = vol2v(v)
% FORMAT [vol, lat_dim] = arrayutils('vol2v', v)
% v         - ND-array of dimensions [nx ny ...] (image or scalar field)
% vol       - 1D-array (i.e., a vector)
% (lat_dim) - [nx ny ...]
%
% Reshape a ND-volume (image or scalar field) into a vector.
    if nargout >= 2
        lat_dim = size(v);
    end
    vol = v(:);
end

function vol = v2field(v, lat_dim)
% FORMAT vol = arrayutils('v2field', v, lat_dim)
% v         - 1D array (i.e., a vector)
% lat_dim   - Dimensions of the lattice ([nx ny ...])
% vol       - ND-array of dimensions [nx ny ... d] (vector field)
%
% Reshape a vector into a ND-volume (vector field). The dimension of the
% vector space is automatically detected from that of the input vector and
% lat_dim.
    s = size(v, 1);
    vec_dim = s / prod(lat_dim);
    vol = reshape(v, [lat_dim vec_dim]);
end

function [vol, lat_dim, vec_dim] = field2v(v)
% FORMAT [vol, lat_dim, vec_dim] = arrayutils('field2v', v)
% v         - ND-array of dimensions [nx ny ... d] (vector field)
% vol       - 1D-array (i.e., a vector)
% (lat_dim) - Dimensions of the lattice ([nx ny ...])
% (vec_dim) - Dimension of the vector space (d)
%
% Reshape a ND-volume (vector field) into a vector.
    if nargout >= 2
        dim = size(v);
        lat_dim = dim(1:(end-1));
    end
    if nargout >= 3
        vec_dim = dim(end);
    end
    vol = v(:);
end

function out = select_slice(in, indices, dim)
% FORMAT out = arrayutils('select_slice', in, indices, dim)
% in        - Input array.
% indices   - Indices to extract in the chosen dimension.
% dim       - Dimension in which slices are extracted.
%
% Select and extract one or several slices from an ND-array.
    if nargin < 3
        dim = -1;
    end
    s = size(in);
    if dim < 1
        dim = length(s);
    end
    S = struct;
    S.type = '()';
    S.subs = cell(length(s),1);
    for j = 1:length(s)
        if j == dim
            S.subs{j} = indices;
        else
            S.subs{j} = ':';
        end
    end
    out = subsref(in, S);
end

function out = assign_slice(in, slice, indices, dim)
% FORMAT out = arrayutils('assign_slice', in, slice, indices, dim)
% in        - Receiving array.
% slice     - Data slice(s) to assign.
% indices   - Indices where to assign the slice in the chosen dimension.
% dim       - Dimension in which slices are assigned.
%
% Assign one or several slices to an ND-array. If there is a shape mismatch
% between the assigned slices and the larger array, the assigned slices are
% reshaped.
    if nargin < 4
        dim = -1;
    end
    s = size(in);
    if dim < 1
        dim = length(s);
    end
    S = struct;
    S.type = '()';
    S.subs = cell(length(s),1);
    for j = 1:length(s)
        if j == dim
            S.subs{j} = indices;
        else
            S.subs{j} = ':';
        end
    end
    slice_dim = size(subsref(in, S));
    out = subsasgn(in, S, reshape(slice, slice_dim));
end

function B = repmat_along(A, d, n)
% FORMAT B = arrayutils('repmat_along', A, dim, n)
% A - fundamental array to repeat
% d - Direction(s) in which to repeat the array
% n - Number(s) of repeats
%
% Repeat an array along a given direction. 
% - If only one direction is provided, the final array size is
%   [size(A,1) ... size(A,d) * n ... size(A,end)].
% - If several directions are provided, the final array size is
%   [size(A,1) ... size(A,d(1)) * n(1) ... 
%    size(A,d(2)) * n(2) ... size(A,end)].
% - If less number of repeats than directions are provided, the last number
%   of repeats is used for the remaining directions.
    s = size(A);
    d = d(:);
    n = n(:);
    n = padarray(n, numel(d) - numel(n), 'replicate', 'post');
    ndim = max([length(s); d]);
    repmat_dim = ones(ndim, 1);
    repmat_dim(d) = n;
    B = repmat(A, repmat_dim');
end