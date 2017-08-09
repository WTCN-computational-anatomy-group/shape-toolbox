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

    if nargout == 0
        obj.hr.dim  = size(h);
        obj.hr(:)   = h(:);
        obj.hr.utd  = true;
    end
end