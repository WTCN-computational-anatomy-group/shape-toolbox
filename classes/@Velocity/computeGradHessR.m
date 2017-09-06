function [g, h] = computeGradHessR(obj)
% FORMAT function [g, h] = obj.computeGradHessR()
%
% Compute Gradient and Hessian of the objective function w.r.t.
% residual field. (obj. f. is the negative log-likelihood)

    if obj.Debug, fprintf('* computeGradHessR\n'); end;

    % --- Add Grad/Hess of the likelihood (matching) term
    [g, h] = obj.gradHessMatchingVel(); 
    g = g * obj.SigmaR;
    h = h * obj.SigmaR^2;
    
    % --- Deal with 2d case
    if size(obj.mu, 3) == 1
        h(:,:,:,[3 5 6]) = 0;
        g(:,:,:,3)       = 0;
    end
    
    % --- Additional regularisation in case H is singular
    % Create Identity tensor field
    r = zeros(size(h), 'like', h);
    r(:,:,:,1:3) = 1;
    % Multiply by the pointwise "maximum of diagonal"
    r = 1e-7 * bsxfun(@times, r, max(h(:,:,:,1:3), [], 4));
    h = h + r;
    clear r

    % --- Deal with 2d case
    if size(obj.mu, 3) == 1
        h(:,:,:,3) = 1;
    end

    if nargout == 0
        obj.hr.dim  = size(h);
        obj.hr(:)   = h(:);
        obj.hr.utd  = true;
    end
end