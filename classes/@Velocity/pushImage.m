function [f, c] = pushImage(obj, ipsi, f)
% FORMAT ([pf, (c)]) = obj.pushImage((ipsi), (f))
% (ipsi) - Inverse transform (warps mu to f). [default: obj.ipsi]
% (f)    - Observed image. [default: obj.f]
% Output:
% (pf)   - Pushed image in template space (if no argout, write to obj.pf)
% (c)    - Pushed voxel count             (if no argout, write to obj.pf)
%
% Push the image to template space

    % --- Check if nothing to do
    if nargin == 1 && obj.checkarray('pf') && obj.checkarray('pvox')
        f = obj.pf;
        c = obj.pvox;
        return
    end
    
    % --- Default arguments
    if nargin < 3
        f = obj.f;
        if nargin < 2
            obj.reconstructIPsi();
            ipsi = obj.ipsi;
        end
    end
    
    % --- Check that all arrays are ready to be used
    if ~obj.checkarray(ipsi) || ~obj.checkarray(f)
        if obj.Debug
            warning('Cannot push image: missing arrays');
        end
        return
    end
    if obj.Debug, fprintf('* pushImage\n'); end;
    
    % --- Load data
    f = single(numeric(f));
    ipsi = single(numeric(ipsi));
        
    % --- Push
    [f, c] = spm_diffeo('pushc', f, ipsi);
    f(~isfinite(f)) = nan;
    c(~isfinite(c)) = nan;
    
    % --- Write or return data
    if nargout == 0
        obj.pf.dim      = size(f);
        obj.pf(:)       = f(:);
        obj.pvox.dim    = size(c);
        obj.pvox(:)     = c(:);
        obj.statusChanged('pf', 'pvox');
    end
end
