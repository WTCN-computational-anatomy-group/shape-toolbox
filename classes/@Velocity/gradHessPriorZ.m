function [g,h] = gradHessPriorZ(obj, z, regz)
% FORMAT [g,h] = obj.gradHessPriorZ((g), (h))
% obj - A velocity object
% (g) - Gradient of the objective function w.r.t. (expected) latent 
%       coordinates. Incremented by this function.
% (h) - Hessian of the objective function w.r.t. (expected) latent 
%       coordinates. Incremented by this function.
%
% Object fields [w, regz] are used by this function 
% and should thus be correctly set before call.
%
% Add the prior part to Gradient and Hessian
% This function should only be called by computeGradHess.

    if nargin < 3
        obj.computeRegZ();
        regz = obj.regz;
        if nargin < 2
            z = obj.z;
        end
    end

    % --- Check all arrays are ready to be used
    if ~obj.checkarray(z) || ~obj.checkarray(regz)
        if obj.Debug
            warning('Cannot compute gradient and hessian: missing arrays');
        end
        return
    end
    if obj.Debug, fprintf('* gradHessPriorZ\n'); end;
    
    % --- Dim info
    dim_latent = size(z, 1);
    
    % --- Load data
    z    = single(numeric(z));
    regz = single(numeric(regz));
    
    % --- Compute grad and hess
    g = zeros([dim_latent 1]);
    g = g + regz * z;
    if nargout > 1
        h = zeros(dim_latent);
        h = h + regz;
    end

end