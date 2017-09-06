function c = pointwise3(a, b, op)
% FORMAT c = pointwise(a, (b), (op))
% a  - [nx nz nz 3 (3)] ND-array of vectors/matrices
% b  - [nx nz nz 3 (3)] ND-array of vectors/matrices
% op - Operation(s) to apply to an element of a before multiplying it 
%      with an element of b:
%      't' (transpose), 'i' (invert), 'd' (determinant).
%      (they can be concatenated: 'itd' means det(inv(a).')
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
    if nargin < 3
        if nargin < 2
            op = '';
            binary = false;
        elseif ischar(b)
            op = b;
            binary = false;
        else
            op = '';
            binary = true;
        end
    else
        binary = true;
    end
    
    % --- Check dimensions
    ismat = struct;
    dim   = struct;
    sym   = struct;
    dim.a = size(a);
    sym.a = size(a, 4) > 3;
    ismat.a = sym.a || size(a, 5) > 1;
    if length(dim.a) < 4
        error('Inputs must be vector or tensor fields')
    end
    if binary
        dim.b = size(b);
        sym.b = size(b, 4) > 3;
        ismat.b = sym.b || size(b, 5) > 1;
        if length(dim.b) < 4
            error('Inputs must be vector or tensor fields')
        end
    end
    if binary && ~issame(dim.a(1:3), dim.b(1:3))
        error('Input fields must have the same lattice')
    end
    
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
                if syma
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
        case 'i'
            if binary && ismat.a && ismat.b
                if sym.a && sym.b
                    c = sasb(isa(a), b);
                elseif sym.a
                    c = samb(isa(a), b);
                elseif sym.b
                    c = masb(ima(a), b);
                else
                    c = mamb(ima(a), b);
                end
            elseif binary && ismat.a && ~ismat.b
                if sym.a
                    c = savb(isa(a), b);
                else
                    c = mavb(ima(a), b);
                end
            elseif ~binary && ismat.a
                if sym.a
                    c = isa(a);
                else
                    c = ima(a);
                end
            else
                error('This case does not make sense')
            end
        case 'd'
            if binary && ismat.a
                if sym.a
                    c = bsxfun(@times, dsa(a), b);
                else
                    c = bsxfun(@times, dma(a), b);
                end
            elseif ~binary && ismat.a
                if sym.a
                    c = dsa(a);
                else
                    c = dma(a);
                end
            else
                error('This case does not make sense')
            end
        case {'it', 'ti'}
            if binary && ismat.a && ismat.b
                if sym.a && sym.b
                    c = sasb(isa(a), b);
                elseif sym.a
                    c = samb(isa(a), b);
                elseif sym.b
                    c = tmasb(ima(a), b);
                else
                    c = tmamb(ima(a), b);
                end
            elseif  binary && ismat.a && ~ismat.b
                if sym.a
                    c = savb(isa(a), b);
                else
                    c = tmavb(ima(a), b);
                end
            elseif ~binary && ismat.a
                if sym.a
                    c = isa(a);
                else
                    c = tma(ima(a));
                end
            else
                error('This case does not make sense')
            end
        case {'itd', 'tid'}
            if binary && ismat.a
                if sym.a
                    c = bsxfun(@times, b, dsa(a));
                else
                    c = bsxfun(@ldivide, b, dma(a));
                end
            elseif ~binary && ismat.a
                if sym.a
                    c = dsa(a);
                else
                    c = 1 ./ dma(a);
                end
            else
                error('This case does not make sense')
            end
        otherwise
            error('The provided operation does not make sense')
    end
end

function c = mamb(a, b)
% a - 3x3 tensor field
% b - 3x3 tensor field
% c - 3x3 tensor field
% Returns pointwise a * b
    dim = size(a);
    c = zeros([dim(1:3) 3 3], 'like', a);
    
    c(:,:,:,1,1) = a(:,:,:,1,1) .* b(:,:,:,1,1) + a(:,:,:,1,2) .* b(:,:,:,2,1) + a(:,:,:,1,3) .* b(:,:,:,3,1);
    c(:,:,:,1,2) = a(:,:,:,1,1) .* b(:,:,:,1,2) + a(:,:,:,1,2) .* b(:,:,:,2,2) + a(:,:,:,1,3) .* b(:,:,:,3,2);
    c(:,:,:,1,3) = a(:,:,:,1,1) .* b(:,:,:,1,3) + a(:,:,:,1,2) .* b(:,:,:,2,3) + a(:,:,:,1,3) .* b(:,:,:,3,3);
    
    c(:,:,:,2,1) = a(:,:,:,2,1) .* b(:,:,:,1,1) + a(:,:,:,2,2) .* b(:,:,:,2,1) + a(:,:,:,2,3) .* b(:,:,:,3,1);
    c(:,:,:,2,2) = a(:,:,:,2,1) .* b(:,:,:,1,2) + a(:,:,:,2,2) .* b(:,:,:,2,2) + a(:,:,:,2,3) .* b(:,:,:,3,2);
    c(:,:,:,2,3) = a(:,:,:,2,1) .* b(:,:,:,1,3) + a(:,:,:,2,2) .* b(:,:,:,2,3) + a(:,:,:,2,3) .* b(:,:,:,3,3);
    
    c(:,:,:,3,1) = a(:,:,:,3,1) .* b(:,:,:,1,1) + a(:,:,:,3,2) .* b(:,:,:,2,1) + a(:,:,:,3,3) .* b(:,:,:,3,1);
    c(:,:,:,3,2) = a(:,:,:,3,1) .* b(:,:,:,1,2) + a(:,:,:,3,2) .* b(:,:,:,2,2) + a(:,:,:,3,3) .* b(:,:,:,3,2);
    c(:,:,:,3,3) = a(:,:,:,3,1) .* b(:,:,:,1,3) + a(:,:,:,3,2) .* b(:,:,:,2,3) + a(:,:,:,3,3) .* b(:,:,:,3,3);
end

function c = masb(a, b)
% a - 3x3 tensor field
% b - 3x3 symmetric tensor field = 6d vector field
% c - 3x3 tensor field
% Returns pointwise a * b
    dim = size(a);
    c = zeros([dim(1:3) 3 3], 'like', a);
    
    ind = symIndices(size(b, 4));
    
    c(:,:,:,1,1) = a(:,:,:,1,1) .* b(:,:,:,ind(1,1)) + a(:,:,:,1,2) .* b(:,:,:,ind(2,1)) + a(:,:,:,1,3) .* b(:,:,:,ind(3,1));
    c(:,:,:,1,2) = a(:,:,:,1,1) .* b(:,:,:,ind(1,2)) + a(:,:,:,1,2) .* b(:,:,:,ind(2,2)) + a(:,:,:,1,3) .* b(:,:,:,ind(3,2));
    c(:,:,:,1,3) = a(:,:,:,1,1) .* b(:,:,:,ind(1,3)) + a(:,:,:,1,2) .* b(:,:,:,ind(2,3)) + a(:,:,:,1,3) .* b(:,:,:,ind(3,3));
    
    c(:,:,:,2,1) = a(:,:,:,2,1) .* b(:,:,:,ind(1,1)) + a(:,:,:,2,2) .* b(:,:,:,ind(2,1)) + a(:,:,:,2,3) .* b(:,:,:,ind(3,1));
    c(:,:,:,2,2) = a(:,:,:,2,1) .* b(:,:,:,ind(1,2)) + a(:,:,:,2,2) .* b(:,:,:,ind(2,2)) + a(:,:,:,2,3) .* b(:,:,:,ind(3,2));
    c(:,:,:,2,3) = a(:,:,:,2,1) .* b(:,:,:,ind(1,3)) + a(:,:,:,2,2) .* b(:,:,:,ind(2,3)) + a(:,:,:,2,3) .* b(:,:,:,ind(3,3));
    
    c(:,:,:,3,1) = a(:,:,:,3,1) .* b(:,:,:,ind(1,1)) + a(:,:,:,3,2) .* b(:,:,:,ind(2,1)) + a(:,:,:,3,3) .* b(:,:,:,ind(3,1));
    c(:,:,:,3,2) = a(:,:,:,3,1) .* b(:,:,:,ind(1,2)) + a(:,:,:,3,2) .* b(:,:,:,ind(2,2)) + a(:,:,:,3,3) .* b(:,:,:,ind(3,2));
    c(:,:,:,3,3) = a(:,:,:,3,1) .* b(:,:,:,ind(1,3)) + a(:,:,:,3,2) .* b(:,:,:,ind(2,3)) + a(:,:,:,3,3) .* b(:,:,:,ind(3,3));
end

function c = samb(a, b)
% a - 3x3 symmetric tensor field = 6d vector field
% b - 3x3 tensor field
% c - 3x3 tensor field
% Returns pointwise a * b
    dim = size(a);
    c = zeros([dim(1:3) 3 3], 'like', a);
    
    ind = symIndices(size(a, 4));
    
    c(:,:,:,1,1) = a(:,:,:,ind(1,1)) .* b(:,:,:,1,1) + a(:,:,:,ind(1,2)) .* b(:,:,:,2,1) + a(:,:,:,ind(1,3)) .* b(:,:,:,3,1);
    c(:,:,:,1,2) = a(:,:,:,ind(1,1)) .* b(:,:,:,1,2) + a(:,:,:,ind(1,2)) .* b(:,:,:,2,2) + a(:,:,:,ind(1,3)) .* b(:,:,:,3,2);
    c(:,:,:,1,3) = a(:,:,:,ind(1,1)) .* b(:,:,:,1,3) + a(:,:,:,ind(1,2)) .* b(:,:,:,2,3) + a(:,:,:,ind(1,3)) .* b(:,:,:,3,3);
    
    c(:,:,:,2,1) = a(:,:,:,ind(2,1)) .* b(:,:,:,1,1) + a(:,:,:,ind(2,2)) .* b(:,:,:,2,1) + a(:,:,:,ind(2,3)) .* b(:,:,:,3,1);
    c(:,:,:,2,2) = a(:,:,:,ind(2,1)) .* b(:,:,:,1,2) + a(:,:,:,ind(2,2)) .* b(:,:,:,2,2) + a(:,:,:,ind(2,3)) .* b(:,:,:,3,2);
    c(:,:,:,2,3) = a(:,:,:,ind(2,1)) .* b(:,:,:,1,3) + a(:,:,:,ind(2,2)) .* b(:,:,:,2,3) + a(:,:,:,ind(2,3)) .* b(:,:,:,3,3);
    
    c(:,:,:,3,1) = a(:,:,:,ind(3,1)) .* b(:,:,:,1,1) + a(:,:,:,ind(3,2)) .* b(:,:,:,2,1) + a(:,:,:,ind(3,3)) .* b(:,:,:,3,1);
    c(:,:,:,3,2) = a(:,:,:,ind(3,1)) .* b(:,:,:,1,2) + a(:,:,:,ind(3,2)) .* b(:,:,:,2,2) + a(:,:,:,ind(3,3)) .* b(:,:,:,3,2);
    c(:,:,:,3,3) = a(:,:,:,ind(3,1)) .* b(:,:,:,1,3) + a(:,:,:,ind(3,2)) .* b(:,:,:,2,3) + a(:,:,:,ind(3,3)) .* b(:,:,:,3,3);
end

function c = sasb(a, b)
% a - 3x3 symmetric tensor field = 6d vector field
% b - 3x3 symmetric tensor field = 6d vector field
% c - 3x3 tensor field
% Returns pointwise a * b
    dim = size(a);
    c = zeros([dim(1:3) 3 3], 'like', a);
    
    ind = symIndices(size(a, 4));
    
    c(:,:,:,1,1) = a(:,:,:,ind(1,1)) .* b(:,:,:,ind(1,1)) + a(:,:,:,ind(1,2)) .* b(:,:,:,ind(2,1)) + a(:,:,:,ind(1,3)) .* b(:,:,:,ind(3,1));
    c(:,:,:,1,2) = a(:,:,:,ind(1,1)) .* b(:,:,:,ind(1,2)) + a(:,:,:,ind(1,2)) .* b(:,:,:,ind(2,2)) + a(:,:,:,ind(1,3)) .* b(:,:,:,ind(3,2));
    c(:,:,:,1,3) = a(:,:,:,ind(1,1)) .* b(:,:,:,ind(1,3)) + a(:,:,:,ind(1,2)) .* b(:,:,:,ind(2,3)) + a(:,:,:,ind(1,3)) .* b(:,:,:,ind(3,3));
    
    c(:,:,:,2,1) = a(:,:,:,ind(2,1)) .* b(:,:,:,ind(1,1)) + a(:,:,:,ind(2,2)) .* b(:,:,:,ind(2,1)) + a(:,:,:,ind(2,3)) .* b(:,:,:,ind(3,1));
    c(:,:,:,2,2) = a(:,:,:,ind(2,1)) .* b(:,:,:,ind(1,2)) + a(:,:,:,ind(2,2)) .* b(:,:,:,ind(2,2)) + a(:,:,:,ind(2,3)) .* b(:,:,:,ind(3,2));
    c(:,:,:,2,3) = a(:,:,:,ind(2,1)) .* b(:,:,:,ind(1,3)) + a(:,:,:,ind(2,2)) .* b(:,:,:,ind(2,3)) + a(:,:,:,ind(2,3)) .* b(:,:,:,ind(3,3));
    
    c(:,:,:,3,1) = a(:,:,:,ind(3,1)) .* b(:,:,:,ind(1,1)) + a(:,:,:,ind(3,2)) .* b(:,:,:,ind(2,1)) + a(:,:,:,ind(3,3)) .* b(:,:,:,ind(3,1));
    c(:,:,:,3,2) = a(:,:,:,ind(3,1)) .* b(:,:,:,ind(1,2)) + a(:,:,:,ind(3,2)) .* b(:,:,:,ind(2,2)) + a(:,:,:,ind(3,3)) .* b(:,:,:,ind(3,2));
    c(:,:,:,3,3) = a(:,:,:,ind(3,1)) .* b(:,:,:,ind(1,3)) + a(:,:,:,ind(3,2)) .* b(:,:,:,ind(2,3)) + a(:,:,:,ind(3,3)) .* b(:,:,:,ind(3,3));
end

function c = mavb(a, b)
% a - 3x3 tensor field
% b - 3d  vector field
% b - 3d  vector field
% Returns pointwise a * b
    dim = size(a);
    c = zeros([dim(1:3) 3], 'like', a);
    
    c(:,:,:,1) = a(:,:,:,1,1) .* b(:,:,:,1) + a(:,:,:,1,2) .* b(:,:,:,2) + a(:,:,:,1,3) .* b(:,:,:,3);
    c(:,:,:,2) = a(:,:,:,2,1) .* b(:,:,:,1) + a(:,:,:,2,2) .* b(:,:,:,2) + a(:,:,:,2,3) .* b(:,:,:,3);
    c(:,:,:,3) = a(:,:,:,3,1) .* b(:,:,:,1) + a(:,:,:,3,2) .* b(:,:,:,2) + a(:,:,:,3,3) .* b(:,:,:,3);
end

function c = savb(a, b)
% a - 3x3 symmetric tensor field = 6d vector field
% b - 3d  vector field
% b - 3d  vector field
% Returns pointwise a * b
    dim = size(a);
    c = zeros([dim(1:3) 3], 'like', a);
    
    ind = symIndices(size(a, 4));
    
    c(:,:,:,1) = a(:,:,:,ind(1,1)) .* b(:,:,:,1) + a(:,:,:,ind(1,2)) .* b(:,:,:,2) + a(:,:,:,ind(1,3)) .* b(:,:,:,3);
    c(:,:,:,2) = a(:,:,:,ind(2,1)) .* b(:,:,:,1) + a(:,:,:,ind(2,2)) .* b(:,:,:,2) + a(:,:,:,ind(2,3)) .* b(:,:,:,3);
    c(:,:,:,3) = a(:,:,:,ind(3,1)) .* b(:,:,:,1) + a(:,:,:,ind(3,2)) .* b(:,:,:,2) + a(:,:,:,ind(3,3)) .* b(:,:,:,3);
end

function c = vamb(a, b)
% a - 3d  vector field
% b - 3x3 tensor field
% b - 3d  vector field
% Returns pointwise a.' * b
    dim = size(a);
    c = zeros([dim(1:3) 3], 'like', a);
    
    c(:,:,:,1) = a(:,:,:,1) .* b(:,:,:,1,1) + a(:,:,:,2) .* b(:,:,:,2,1) + a(:,:,:,3) .* b(:,:,:,3,1);
    c(:,:,:,2) = a(:,:,:,1) .* b(:,:,:,1,2) + a(:,:,:,2) .* b(:,:,:,2,2) + a(:,:,:,3) .* b(:,:,:,3,2);
    c(:,:,:,3) = a(:,:,:,1) .* b(:,:,:,1,3) + a(:,:,:,2) .* b(:,:,:,2,3) + a(:,:,:,3) .* b(:,:,:,3,3);
end

function c = vasb(a, b)
% a - 3d  vector field
% b - 3x3 symmetric tensor field = 6d vector field
% b - 3d  vector field
% Returns pointwise a.' * b
    dim = size(a);
    c = zeros([dim(1:3) 3], 'like', a);
    
    ind = symIndices(size(b, 4));
    
    c(:,:,:,1) = a(:,:,:,1) .* b(:,:,:,ind(1,1)) + a(:,:,:,2) .* b(:,:,:,ind(2,1)) + a(:,:,:,3) .* b(:,:,:,ind(3,1));
    c(:,:,:,2) = a(:,:,:,1) .* b(:,:,:,ind(1,2)) + a(:,:,:,2) .* b(:,:,:,ind(2,2)) + a(:,:,:,3) .* b(:,:,:,ind(3,2));
    c(:,:,:,3) = a(:,:,:,1) .* b(:,:,:,ind(1,3)) + a(:,:,:,2) .* b(:,:,:,ind(2,3)) + a(:,:,:,3) .* b(:,:,:,ind(3,3));
end

function c = vavb(a, b)
% a - 3d vector field
% b - 3d vector field
% b - 3d vector field
% Returns pointwise a.' * b
    c = a(:,:,:,1) .* b(:,:,:,1) + a(:,:,:,2) .* b(:,:,:,2) + a(:,:,:,3) .* b(:,:,:,3);
end

function c = tmamb(a, b)
% a - 3x3 tensor field
% b - 3x3 tensor field
% c - 3x3 tensor field
% Returns pointwise a.' * b
    dim = size(a);
    c = zeros([dim(1:3) 3 3], 'like', a);
    
    c(:,:,:,1,1) = a(:,:,:,1,1) .* b(:,:,:,1,1) + a(:,:,:,2,1) .* b(:,:,:,2,1) + a(:,:,:,3,1) .* b(:,:,:,3,1);
    c(:,:,:,1,2) = a(:,:,:,1,1) .* b(:,:,:,1,2) + a(:,:,:,2,1) .* b(:,:,:,2,2) + a(:,:,:,3,1) .* b(:,:,:,3,2);
    c(:,:,:,1,3) = a(:,:,:,1,1) .* b(:,:,:,1,3) + a(:,:,:,2,1) .* b(:,:,:,2,3) + a(:,:,:,3,1) .* b(:,:,:,3,3);
    
    c(:,:,:,2,1) = a(:,:,:,1,2) .* b(:,:,:,1,1) + a(:,:,:,2,2) .* b(:,:,:,2,1) + a(:,:,:,3,2) .* b(:,:,:,3,1);
    c(:,:,:,2,2) = a(:,:,:,1,2) .* b(:,:,:,1,2) + a(:,:,:,2,2) .* b(:,:,:,2,2) + a(:,:,:,3,2) .* b(:,:,:,3,2);
    c(:,:,:,2,3) = a(:,:,:,1,2) .* b(:,:,:,1,3) + a(:,:,:,2,2) .* b(:,:,:,2,3) + a(:,:,:,3,2) .* b(:,:,:,3,3);
    
    c(:,:,:,3,1) = a(:,:,:,1,3) .* b(:,:,:,1,1) + a(:,:,:,2,3) .* b(:,:,:,2,1) + a(:,:,:,3,3) .* b(:,:,:,3,1);
    c(:,:,:,3,2) = a(:,:,:,1,3) .* b(:,:,:,1,2) + a(:,:,:,2,3) .* b(:,:,:,2,2) + a(:,:,:,3,3) .* b(:,:,:,3,2);
    c(:,:,:,3,3) = a(:,:,:,1,3) .* b(:,:,:,1,3) + a(:,:,:,2,3) .* b(:,:,:,2,3) + a(:,:,:,3,3) .* b(:,:,:,3,3);
end

function c = tma(a)
% a - 3x3 tensor field
% c - 3x3 tensor field
% Returns pointwise a.'
    dim = size(a);
    c = zeros([dim(1:3) 3 3], 'like', a);
    
    c(:,:,:,1,1) = a(:,:,:,1,1);
    c(:,:,:,1,2) = a(:,:,:,2,1);
    c(:,:,:,1,3) = a(:,:,:,3,1);
    
    c(:,:,:,2,1) = a(:,:,:,1,2);
    c(:,:,:,2,2) = a(:,:,:,2,2);
    c(:,:,:,2,3) = a(:,:,:,3,2);
    
    c(:,:,:,3,1) = a(:,:,:,1,3);
    c(:,:,:,3,2) = a(:,:,:,2,3);
    c(:,:,:,3,3) = a(:,:,:,3,3);
end

function c = ima(a)
% a - 3x3 tensor field
% c - 3x3 tensor field
% Returns pointwise inv(a)
    dim = size(a);
    c = zeros([dim(1:3) 3 3], 'like', a);
    
    % - Cofactors
    
    c(:,:,:,1,1) = a(:,:,:,2,2) .* a(:,:,:,3,3) - a(:,:,:,2,3) .* a(:,:,:,3,2);
    c(:,:,:,1,2) = a(:,:,:,1,3) .* a(:,:,:,3,2) - a(:,:,:,1,2) .* a(:,:,:,3,3);
    c(:,:,:,1,3) = a(:,:,:,1,2) .* a(:,:,:,2,3) - a(:,:,:,1,3) .* a(:,:,:,2,2);
    
    c(:,:,:,2,1) = a(:,:,:,2,3) .* a(:,:,:,3,1) - a(:,:,:,2,1) .* a(:,:,:,3,3);
    c(:,:,:,2,2) = a(:,:,:,1,1) .* a(:,:,:,3,3) - a(:,:,:,1,3) .* a(:,:,:,3,1);
    c(:,:,:,2,3) = a(:,:,:,1,3) .* a(:,:,:,2,1) - a(:,:,:,1,1) .* a(:,:,:,2,3);
    
    c(:,:,:,3,1) = a(:,:,:,2,1) .* a(:,:,:,3,2) - a(:,:,:,2,2) .* a(:,:,:,3,1);
    c(:,:,:,3,2) = a(:,:,:,1,2) .* a(:,:,:,3,1) - a(:,:,:,1,1) .* a(:,:,:,3,2);
    c(:,:,:,3,3) = a(:,:,:,1,1) .* a(:,:,:,2,2) - a(:,:,:,1,2) .* a(:,:,:,2,1);
    
    % - Determinant
    
    c = bsxfun(@rdivide, c, a(:,:,:,1,1) .* c(:,:,:,1,1) + ...
                             a(:,:,:,1,2) .* c(:,:,:,2,1) + ...
                             a(:,:,:,1,3) .* c(:,:,:,3,1) );
end

function c = isa(a)
% a - 3x3 symmetric tensor field = 6d vector field
% c - 3x3 symmetric tensor field = 6d vector field
% Returns pointwise inv(a)
    dim = size(a);
    c = zeros([dim(1:3) 6], 'like', a);
    
    ind = symIndices(size(a, 4));
    
    % - Cofactors
    
    c(:,:,:,ind(1,1)) = a(:,:,:,ind(2,2)) .* a(:,:,:,ind(3,3)) - a(:,:,:,ind(2,3)) .* a(:,:,:,ind(3,2));
    c(:,:,:,ind(1,2)) = a(:,:,:,ind(1,3)) .* a(:,:,:,ind(3,2)) - a(:,:,:,ind(1,2)) .* a(:,:,:,ind(3,3));
    c(:,:,:,ind(1,3)) = a(:,:,:,ind(1,2)) .* a(:,:,:,ind(2,3)) - a(:,:,:,ind(1,3)) .* a(:,:,:,ind(2,2));
    
    c(:,:,:,ind(2,2)) = a(:,:,:,ind(1,1)) .* a(:,:,:,ind(3,3)) - a(:,:,:,ind(1,3)) .* a(:,:,:,ind(3,1));
    
    c(:,:,:,ind(3,2)) = a(:,:,:,ind(1,2)) .* a(:,:,:,ind(3,1)) - a(:,:,:,ind(1,1)) .* a(:,:,:,ind(3,2));
    
    % - Determinant
    
    c = bsxfun(@rdivide, c, a(:,:,:,ind(1,1)) .* c(:,:,:,ind(1,1)) + ...
                            a(:,:,:,ind(1,2)) .* c(:,:,:,ind(2,1)) + ...
                            a(:,:,:,ind(1,3)) .* c(:,:,:,ind(3,1)) );
end

function c = dma(a)
% a - 3x3 tensor field
% c - scalar field
% Returns pointwise det(a)
    c = a(:,:,:,1,1) .* ( a(:,:,:,2,2) .* a(:,:,:,3,3) - a(:,:,:,2,3) .* a(:,:,:,3,2) ) + ...
        a(:,:,:,1,2) .* ( a(:,:,:,2,3) .* a(:,:,:,3,1) - a(:,:,:,2,1) .* a(:,:,:,3,3) ) + ...
        a(:,:,:,1,3) .* ( a(:,:,:,2,1) .* a(:,:,:,3,2) - a(:,:,:,2,2) .* a(:,:,:,3,1) );
end


function c = dsa(a)
% a - 3x3 symmetric tensor field = 6d vector field
% c - scalar field
% Returns pointwise det(a)
    ind = symIndices(size(a, 4));
    c = a(:,:,:,ind(1,1)) .* ( a(:,:,:,ind(2,2)) .* a(:,:,:,ind(3,3)) - a(:,:,:,ind(2,3)) .* a(:,:,:,ind(3,2)) ) + ...
        a(:,:,:,ind(1,2)) .* ( a(:,:,:,ind(2,3)) .* a(:,:,:,ind(3,1)) - a(:,:,:,ind(2,1)) .* a(:,:,:,ind(3,3)) ) + ...
        a(:,:,:,ind(1,3)) .* ( a(:,:,:,ind(2,1)) .* a(:,:,:,ind(3,2)) - a(:,:,:,ind(2,2)) .* a(:,:,:,ind(3,1)) );
end