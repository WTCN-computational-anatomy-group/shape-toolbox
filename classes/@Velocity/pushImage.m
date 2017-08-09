function [f, c] = pushImage(obj, iphi, f)
% FORMAT ([pf, (c)]) = obj.pushImage((iphi), (f))
% (iphi) - Inverse transform (warps mu to f). [default: obj.iphi]
% (f)    - Observed image. [default: obj.f]
% Output:
% (pf)   - Pushed image in template space (if no argout, write to obj.pf)
% (c)    - Pushed voxel count             (if no argout, write to obj.pf)
%
% Push the image to template space

    % --- Check if nothing to do
    if nargout == 0 && obj.checkarray('pf') && obj.checkarray('pvox')
        return
    end
    
    % --- Default arguments
    if nargin < 3
        f = obj.f;
        if nargin < 2
            obj.exponentiateVelocity('iphi');
            iphi = obj.iphi;
        end
    end
    
    % --- Check that all arrays are ready to be used
    if ~obj.checkarray(iphi) || ~obj.checkarray(f)
        if obj.Debug
            warning('Cannot push image: missing arrays');
        end
        return
    end
    if obj.Debug, fprintf('* pushImage\n'); end;
    
    % --- Load data
    if isa(f, 'file_array')
        f = single(numeric(f));
    end
    if isa(iphi, 'file_array')
        iphi = single(numeric(iphi));
    end
        
    % --- Push
    [f, c] = spm_diffeo('pushc', f, iphi);
    f(~isfinite(f)) = nan;
    c(~isfinite(c)) = nan;
    
    % --- Write or return data
    if nargout == 0
        obj.pf.dim      = size(f);
        obj.pf(:)       = f(:);
        obj.pvox.dim    = size(c);
        obj.pvox(:)     = c(:);
        obj.statusChanged('pf', 'pvox');
        obj.utd.pf      = true;
        obj.utd.pvox    = true;
    end
end
