function c = pointwise(a, varargin)
% FORMAT c = pointwise(a, (b), (op), (sym))
% a  - [nx nz nz A (B)] ND-array of vectors/matrices
% b  - [nx nz nz C (D)] ND-array of vectors/matrices
% op - Operation(s) to apply to an element of a before multiplying it 
%      with an element of b:
%      't' (transpose).
% sym - If true, at least one of a and b is symmetric.
%
% Performs matrix operations at each point of a tensor field.
% If b is provided, at each point:
%    c = op(a) * b
% Else, at each point:
%    c = op(a)
%
% Note that
% - if one (or both) entry is a vector, the adequate transpose
%   operation is usually automatically detected.
% - a and b can be fields of sparse square matrices. In this case, they
%   have the form of a 6 element vector, which is automatically detected.

    % --- Default arguments
    b = [];
    issym = false;
    op = '';
    while ~isempty(varargin)
        if ischar(varargin{1})
            op = varargin{1};
        elseif isscalar(varargin{1})
            issym = varargin{1};
        else
            b = varargin{1};
        end
        varargin = varargin(2:end);
    end
    binary = ~isempty(b);
    
    % --- Check dimensions
    ismat = struct;
    dim   = struct;
    sym   = struct;
    dim.a = [size(a) 1 1];
    if ~binary
        sym.a = issym;
    end
    if binary
        dim.b = [size(b) 1 1];
        if issym
            if size(a,5) == 1 && size(b,5) > 1
                sym.a = true;
                sym.b = false;
            elseif size(a,5) > 1 && size(b,5) == 1
                sym.b = true;
                sym.a = false;
            else
                sym.a = false;
                sym.b = false;
                if size(a,4) >=  size(b,4)
                    sym.a = true;
                end
                if size(b,4) >= size(a,4)
                    sym.b = true;
                end
            end
        else
            sym.a = false;
            sym.b = false;
        end
    end
    if binary && ~issame(dim.a(1:3), dim.b(1:3))
        error('Input fields must have the same lattice')
    end
    ismat.a = sym.a || size(a, 5) > 1;
    ismat.b = sym.b || size(b, 5) > 1;
    
    % --- Deal with cases
    switch lower(op)
        case ''
            if binary && ismat.a && ismat.b
                if sym.a && sym.b
                    c = sasb(a, b);
                elseif sym.a
                    c = samb(a, b);
                elseif sym.b
                    c = masb(a, b);
                else
                    c = mamb(a, b);
                end
            elseif binary && ismat.a && ~ismat.b
                if sym.a
                    c = savb(a, b);
                else
                    c = mavb(a, b);
                end
            elseif binary && ~ismat.a && ismat.b
                if sym.b
                    c = vasb(a, b);
                else
                    c = vamb(a, b);
                end
            elseif  binary && ~ismat.a && ~ismat.b
                c = vavb(a, b);
            else
                c = a;
            end
        case 't'
            if binary && ismat.a && ismat.b
                if sym.a && sym.b
                    c = sasb(a, b);
                elseif sym.b
                    c = masb(a, b);
                elseif sym.a
                    c = samb(a, b);
                else
                    c = tmamb(a, b);
                end
            elseif binary && ismat.a && ~ismat.b
                if sym.a
                    c = savb(a, b);
                else
                    c = tmavb(a, b);
                end
            elseif binary && ~ismat.a && ismat.b
                if sym.b
                    c = vasb(a, b);
                else
                    c = vamb(a, b);
                end
            elseif  binary && ~ismat.a && ~ismat.b
                c = vavb(a, b);
            elseif ismat.a
                if sym.a
                    c = a;
                else
                    c = tma(a);
                end
            else
                c = a;
            end
        otherwise
            error('The provided operation is not supported')
    end
end

function c = mamb(a, b)
% a - 3x3 tensor field
% b - 3x3 tensor field
% c - 3x3 tensor field
% Returns pointwise a * b

    if size(a,5) ~= size(b,4)
        error('Matrix dimensions not consistant: %d and %d', size(a,5), size(b,4))
    end
    
    dim = [size(a) 1 1];
    c = zeros([dim(1:3) size(a,4) size(b,5)], 'like', a);
    
    for i=1:size(a,4)
        for j=1:size(b,5)
            for k=1:size(a,5)
                c(:,:,:,i,j) = c(:,:,:,i,j) + a(:,:,:,i,k) .* b(:,:,:,k,j);
            end
        end
    end
end

function c = masb(a, b)
% a - 3x3 tensor field
% b - 3x3 symmetric tensor field = 6d vector field
% c - 3x3 tensor field
% Returns pointwise a * b

    [ind, n] = symIndices(size(b, 4));
    if size(a,5) ~= n
        error('Matrix dimensions not consistant: %d and %d', size(a,5), n)
    end
    
    dim = [size(a) 1 1];
    c = zeros([dim(1:3) size(a,4) n], 'like', a);
    
    for i=1:size(a,4)
        for j=1:n
            for k=1:n
                c(:,:,:,i,j) = c(:,:,:,i,j) + a(:,:,:,i,k) .* b(:,:,:,ind(k,j));
            end
        end
    end
end

function c = samb(a, b)
% a - 3x3 symmetric tensor field = 6d vector field
% b - 3x3 tensor field
% c - 3x3 tensor field
% Returns pointwise a * b

    ind = symIndices(size(a, 4));
    if n ~= size(b,4)
        error('Matrix dimensions not consistant: %d and %d', n, size(b,4))
    end
    
    dim = [size(a) 1 1];
    c = zeros([dim(1:3) n size(b,5)], 'like', a);
    
    for i=1:n
        for j=1:size(b,5)
            for k=1:n
                c(:,:,:,i,j) = a(:,:,:,ind(i,k)) .* b(:,:,:,k,j);
            end
        end
    end
end

function c = sasb(a, b)
% a - 3x3 symmetric tensor field = 6d vector field
% b - 3x3 symmetric tensor field = 6d vector field
% c - 3x3 tensor field
% Returns pointwise a * b

    
    [inda, na] = symIndices(size(a, 4));
    [indb, nb] = symIndices(size(b, 4));
    if na ~= nb
        error('Matrix dimensions not consistant: %d and %d', na, nb)
    end
    
    dim = [size(a) 1 1];
    c = zeros([dim(1:3) na na], 'like', a);
    
    for i=1:na
        for j=1:na
            for k=1:na
                c(:,:,:,i,j) = c(:,:,:,i,j) + a(:,:,:,inda(i,k)) .* b(:,:,:,indb(k,j));
            end
        end
    end
end

function c = mavb(a, b)
% a - 3x3 tensor field
% b - 3d  vector field
% c - 3d  vector field
% Returns pointwise a * b

    transpose = false;
    n = size(a,4);
    if size(a,5) ~= size(b,4)
        if size(a,4) ~= size(b,4)
            error('Matrix dimensions not consistant: %d and %d', size(a,5), size(b,4))
        else
            n = size(a,5);
            transpose = true;
        end
    end
    
    dim = [size(a) 1 1];
    c = zeros([dim(1:3) n], 'like', a);
    
    if ~transpose
        for i=1:size(a,4)
            for k=1:size(b,4)
                c(:,:,:,i) = c(:,:,:,i) + a(:,:,:,i,k) .* b(:,:,:,k);
            end
        end
    else
        for i=1:size(a,5)
            for k=1:size(b,4)
                c(:,:,:,i) = c(:,:,:,i) + a(:,:,:,k,i) .* b(:,:,:,k);
            end
        end
    end
end

function c = savb(a, b)
% a - 3x3 symmetric tensor field = 6d vector field
% b - 3d  vector field
% c - 3d  vector field
% Returns pointwise a * b

    [ind, n] = symIndices(size(a, 4));
    if n ~= size(b,4)
        error('Matrix dimensions not consistant: %d and %d', n, size(b,4))
    end
    
    dim = [size(a) 1 1];
    c = zeros([dim(1:3) n], 'like', a);
    
    for i=1:n
        for k=1:n
            c(:,:,:,i) = c(:,:,:,i) + a(:,:,:,ind(i,k)) .* b(:,:,:,k);
        end
    end
end

function c = vamb(a, b)
% a - 3d  vector field
% b - 3x3 tensor field
% b - 3d  vector field
% Returns pointwise a.' * b

    transpose = false;
    n = size(b,5);
    if size(a,4) ~= size(b,4)
        if size(a,4) ~= size(b,5)
            error('Matrix dimensions not consistant: %d and %d', size(a,4), size(b,4))
        else
            n = size(b,4);
            transpose = true;
        end
    end
    
    dim = [size(a) 1 1];
    c = zeros([dim(1:3) n], 'like', a);
    
    if ~transpose
        for i=1:size(b,5)
            for k=1:size(a,4)
                c(:,:,:,i) = c(:,:,:,i) + a(:,:,:,k) .* b(:,:,:,k,i);
            end
        end
    else
        for i=1:size(b,4)
            for k=1:size(a,4)
                c(:,:,:,i) = c(:,:,:,i) + a(:,:,:,k) .* b(:,:,:,i,k);
            end
        end
    end
end

function c = vasb(a, b)
% a - 3d  vector field
% b - 3x3 symmetric tensor field = 6d vector field
% b - 3d  vector field
% Returns pointwise a.' * b

    [ind, n] = symIndices(size(b, 4));
    if size(a,4) ~= n
        error('Matrix dimensions not consistant: %d and %d', size(a,4), n)
    end
    
    dim = [size(a) 1 1];
    c = zeros([dim(1:3) n], 'like', a);
    
    for i=1:n
        for k=1:n
            c(:,:,:,i) = c(:,:,:,i) + a(:,:,:,k) .* b(:,:,:,ind(k,i));
        end
    end
end

function c = vavb(a, b)
% a - 3d vector field
% b - 3d vector field
% b - 3d vector field
% Returns pointwise a.' * b
    
    if size(a,4) ~= size(b,4)
        error('Matrix dimensions not consistant: %d and %d', size(a,4), size(b,4))
    end
    
    dim = [size(a) 1 1];
    c = zeros(dim(1:3), 'like', a);
    
    for k=1:size(a,4)
        c = c + a(:,:,:,k) .* b(:,:,:,k);
    end
end

function c = tmamb(a, b)
% a - 3x3 tensor field
% b - 3x3 tensor field
% c - 3x3 tensor field
% Returns pointwise a.' * b

    if size(a,4) ~= size(b,4)
        error('Matrix dimensions not consistant: %d and %d', size(a,4), size(b,4))
    end
    
    dim = [size(a) 1 1];
    c = zeros([dim(1:3) size(a,5) size(b,5)], 'like', a);
    
    for i=1:size(a,5)
        for j=size(b,5)
            for k=1:size(a,4)
                c(:,:,:,i,j) = c(:,:,:,i,j) + a(:,:,:,k,i) .* b(:,:,:,k,j);
            end
        end
    end
end

function c = tmavb(a, b)
% a - 3x3 tensor field
% b - 3d  vector field
% c - 3d  vector field
% Returns pointwise a.' * b
    
    transpose = false;
    n = size(a,5);
    if size(a,4) ~= size(b,4)
        if size(a,5) ~= size(b,4)
            error('Matrix dimensions not consistant: %d and %d', size(a,4), size(b,4))
        else
            n = size(a,4);
            transpose = true;
        end
    end
    
    dim = [size(a) 1 1];
    c = zeros([dim(1:3) n], 'like', numeric(a(1)));
    
    if ~transpose
        for i=1:size(a,5)
            for k=1:size(b,4)
                c(:,:,:,i) = c(:,:,:,i) + a(:,:,:,k,i) .* b(:,:,:,k);
            end
        end
    else
        for i=1:size(a,4)
            for k=1:size(b,4)
                c(:,:,:,i) = c(:,:,:,i) + a(:,:,:,i,k) .* b(:,:,:,k);
            end
        end
    end
end

function c = tma(a)
% a - 3x3 tensor field
% c - 3x3 tensor field
% Returns pointwise a.'

    dim = [size(a) 1 1];
    c = zeros([dim(1:3) size(a,5) size(a,4)], 'like', a);
    
    for i=1:size(a,5)
        for j=1:size(a,4)
            c(:,:,:,i,j) = a(:,:,:,j,i);
        end
    end
end
