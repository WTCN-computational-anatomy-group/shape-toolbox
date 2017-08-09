function [g, h] = computeGradHessZ(obj)
% FORMAT function [g, h] = obj.computeGradHessZ()
%
% Compute Gradient and Hessian of the objective function w.r.t. latent
% coordinates. (obj. f. is the negative log-likelihood)

    if obj.Debug, fprintf('* computeGradHessZ\n'); end;
    
    % --- Initialize Grad/Hess
    dim_latent = size(obj.w, 5);
    g = zeros([dim_latent 1]);
    h = zeros(dim_latent);

    % --- Add Grad/Hess of the likelihood (matching) term
    [g1, h1] = obj.gradHessMatchingZ();
    g = g + g1;
    h = h + h1;
    clear g1 h1

    % --- Add Grad/Hess of the prior term
    [g1, h1] = obj.gradHessPriorZ();
    g = g + g1;
    h = h + h1;
    clear g1 h1
    
    % --- Additional regularisation in case H is singular
    r = 1e-7 * max(diag(h)) * eye(dim_latent);
    h = h + r;
    clear r
    
    if nargout == 0
        obj.hz.dim  = size(h);
        obj.hz(:)   = h(:);
        h           = obj.hz;
        obj.hz.utd  = true;
    end
    
end