function out = defToColor(d, nrm)
% FORMAT c = defToColor(d)
% d   - categorical (4D) image.
%
% Generate an RGB volume from a displacement (e.g. initial velocity)
% volume.

    d = numeric(d);
    if nargin < 2 || ~isfinite(nrm)
        nrm = max(max(max(sum(d.^2, 4)))); % normalising constant
    end

    % Create HSV
    dim = [size(d) 1 1];
    h = zeros(dim(1:3));
    s = sum(d.^2, 4) / nrm;
    for k=1:size(d, 4)
        dk = d(:,:,:,k);
        neg = dk;
        neg(dk > 0) = 0;
        pos = dk;
        pos(dk < 0) = 0;
        h = h + (k-1)*180/3 * abs(pos) + ((k-1)*180/3+180) * abs(neg);
    end
    clear dk neg pos
    h = h./sum(abs(d), 4);
    clear d
    h = mod(h, 360);
    
    out = ones([dim(1:3) 3]);
    out(:,:,:,1) = h/360;
    out(:,:,:,2) = s;
    out = reshape(out, [], dim(3), 3);
    out = hsv2rgb(out);
    out = reshape(out, [dim(1:3) 3]);
                   
end