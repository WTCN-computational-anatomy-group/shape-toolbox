function wf = warpImage(obj, psi, f)
% FORMAT (wf) = obj.warpImage((psi), (f))
% (psi) - Direct transform (warps f to mu). [default: obj.psi]
% (f)   - Image. [default: obj.f]
% (wf)  - Warped image [if no argout: write to obj.wf]
% obj   - Velocity object. This function will use property 'Interpolation'.
%
% Warp the image in template space

    % --- Check if nothing to do
    if nargin == 1 && obj.checkarray('wf')
        return
    end
    
    % --- Default arguments
    if nargin < 3
        f = obj.f;
        if nargin < 2
            obj.reconstructPsi();
            psi = obj.psi;
        end
    end
    
    % --- If all empty, get out
    % I have to check this way because for file_arrays, numel and
    % isempty do not work the same way as for arrays
    if ~obj.checkarray(psi) || ~obj.checkarray(f)
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
    dim = [size(psi) 1 1 1 1];
    dim = dim(1:4);
    lat = dim(1:3);             % Dimension of the output lattice
    
    % --- Prepare output and interpolate
    if nargout == 0
        obj.wf.dim  = [lat nc];
        wf          = obj.wf;
        warp(psi, f, obj.Interpolation, wf);
        obj.statusChanged('wf');
    else
        wf = warp(psi, f, obj.Interpolation);
    end
end