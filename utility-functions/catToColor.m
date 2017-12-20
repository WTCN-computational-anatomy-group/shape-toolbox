function c = catToColor(f, pal)
% FORMAT c = catToColor(f, pal)
% f   - categorical (4D) image.
% pal - palette (Mx3 array or handle to palette function) [hsv]
%
% Generate an RGB volume from a categorical (e.g. responsibilities) volume.

    if nargin < 2
        pal = @hsv;
    end
    if isa(pal, 'function_handle')
        pal = pal(size(f,4));
    end
    
    dim = [size(f) 1 1];
    c   = zeros([dim(1:3) 3]); % output RGB image
    s   = zeros(dim(1:3));     % normalising term
    
    for k=1:dim(4)
        s = s + f(:,:,:,k);
        color = reshape(pal(k,:), [1 1 1 3]);
        c = c + bsxfun(@times, f(:,:,:,k), color);
    end

end