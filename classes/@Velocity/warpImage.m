function wf = warpImage(obj, phi, f)
% FORMAT (wf) = obj.warpImage((phi), (f))
% (phi) - Direct transform (warps f to mu). [default: obj.phi]
% (f)   - Image. [default: obj.f]
% (wf)  - Warped image [if no argout: write to obj.wf]
% obj   - Velocity object. This function will use property 'Interpolation'.
%
% Warp the image in template space

    % --- Check if nothing to do
    if nargout == 0 && obj.checkarray('wf')
        return
    end
    
    % --- Default arguments
    if nargin < 3
        f = obj.f;
        if nargin < 2
            if ~obj.checkarray('phi'), obj.exponentiateVelocity('phi'); end
            phi = obj.phi;
        end
    end
    
    % --- If all empty, get out
    % I have to check this way because for file_arrays, numel and
    % isempty do not work the same way as for arrays
    if ~obj.checkarray(phi) || ~obj.checkarray(f)
        wf = single([]);
        if obj.Debug
            warning('Cannot write image: missing arrays\n')
        end
        return
    end
    if obj.Debug, fprintf('* warpImage\n'); end;
    
    % --- Dim info
    dim = [size(f) 1 1 1 1];
    dim = dim(1:4);
    nc  = dim(3);               % Number of classes/modalities
    dim = [size(phi) 1 1 1 1];
    dim = dim(1:4);
    lat = dim(1:3);             % Dimension of the output lattice
    
    % --- Prepare output and interpolate
    if nargout == 0
        obj.wf.dim  = [lat nc];
        wf          = obj.wf;
        warp(phi, f, obj.Interpolation, wf);
        % --- Up-to-date
        obj.statusChanged('wf');
        obj.utd.wf = true;
    else
        wf = warp(phi, f, obj.Interpolation);
    end
end