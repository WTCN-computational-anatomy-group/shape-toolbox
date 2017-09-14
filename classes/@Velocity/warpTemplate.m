function wa = warpTemplate(obj, ipsi, a)
% FORMAT (wa) = obj.warpTemplate((iphi), (a))
% (ipsi) - Inverse transform (warps mu to f). [default: obj.ipsi]
% (a)    - Template. [default: obj.mu/obj.a]
% (wa)   - Warped template [if no argout: write to obj.wmu/obj.wa]
% obj    - Velocity object. This function uses property 'Interpolation'
%
% Warp the template in image space

    % --- Check if nothing to do
    if nargin == 1 && obj.checkarray('wa')
        return
    end
    
    % --- Default arguments
    if nargin < 3
        a = obj.a;  % There's a direction toward mu for normal/laplace
        if nargin < 2
            obj.reconstructIPsi();
            ipsi = obj.ipsi;
        end
    end
    
    % --- If all empty, get out
    % I have to check this way because for file_arrays, numel and
    % isempty do not work the same way as for arrays
    if ~obj.checkarray(ipsi) || ~obj.checkarray(a)
        if obj.Debug
            warning('Cannot warp template: missing arrays\n')
        end
        wa = single([]);
        return
    end
    if obj.Debug, fprintf('* warpTemplate\n'); end;
    
    % --- Dim info
    dim = [size(a) 1 1 1 1];
    dim = dim(1:4);
    nc  = dim(4);               % Number of classes/modalities
    dim = [size(ipsi) 1 1 1 1];
    dim = dim(1:4);
    lat = dim(1:3);             % Dimension of the output lattice
    
    % --- Prepare output and interpolate
    if nargout == 0
        obj.wa.dim = [lat nc];
        wa         = obj.wa;
        warp(ipsi, a, obj.Interpolation, wa);
        obj.statusChanged('wa');
    else
        wa = warp(ipsi, a, obj.Interpolation);
    end
end