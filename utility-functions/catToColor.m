function c = catToColor(f, pal)
% FORMAT c = catToColor(f, pal)
% f   - categorical (4D) image.
% pal - palette (Mx3 array or handle to palette function) [hsv]
%
% Generate an RGB volume from a categorical (e.g. responsibilities) volume.

    if nargin < 2
        pal = @hsv;
    end
    
    tri = false;
    if numel(size(f)) == 4 && size(f, 3) == 1
        tri = true;
        dim = [size(f) 1 1];
        f = reshape(f, [dim(1:2) dim(4)]);
    end
    if isa(pal, 'function_handle')
        pal = pal(size(f,3));
    end
    
    dim = [size(f) 1 1];
    c   = zeros([dim(1:2) 3]); % output RGB image
    s   = zeros(dim(1:2));     % normalising term
    
    for k=1:dim(3)
        s = s + f(:,:,k);
        color = reshape(pal(k,:), [1 1 3]);
        c = c + bsxfun(@times, f(:,:,k), color);
    end
    c = bsxfun(@rdivide, c, s);

    if tri
        c = reshape(c, [size(c, 1) size(c, 2) 1 size(c, 3)]);
    end
end