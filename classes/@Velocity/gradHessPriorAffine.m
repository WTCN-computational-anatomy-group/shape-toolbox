function [g,h] = gradHessPriorAffine(obj, q, regq)
% FORMAT [g,h] = obj.gradHessPriorAffine((g), (h))
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
        regq = obj.RegAffine;
        if nargin < 2
            q = obj.q;
        end
    end
    
    % --- Check that all arrays are ready to be used
    if ~obj.checkarray(q) || ~obj.checkarray(regq)
        if obj.Debug
            warning('Cannot compute gradient and hessian: missing arrays');
        end
        return
    end
    if obj.Debug, fprintf('* gradHessPriorAffine\n'); end;
    
    
    % --- Dim info
    nq = size(q, 1);
    npq = size(regq, 1);
    
    % --- Load data
    q    = numeric(q);
    regq = numeric(regq);
    % The first basis (rigid) have no prior
    sub  = (nq-npq+1):nq; 
    
    % --- Compute grad and hess
    g = zeros([nq 1]);
    g(sub) = regq * q(sub);
    if nargout > 1
        h = zeros(nq);
        h(sub,sub) = regq;
    end
    
end