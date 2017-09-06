function [g,h] = gradHessPriorZ(obj, z, regz)
% FORMAT [g,h] = obj.gradHessPriorZ((z), (regz))
% obj  - A velocity object
% z    - Latent coordinates
% regz - Inverse covariance matrix of Z (W'LW)
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
    g = regz * z;
    if nargout > 1
        h = regz;
    end

end