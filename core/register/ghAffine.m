% TODO: adapt

function [g, h] = computeGradHessAffine(obj)
% FORMAT function [g, h] = computeGradHessAffine()
%
% Compute Gradient and Hessian of the objective function w.r.t. affine
% parameters. (objective function is the negative log-likelihood).

    if obj.Debug, fprintf('* computeGradHessAffine\n'); end;
    
    % --- Initialize Grad/Hess
    nq = size(obj.q, 1);
    g = zeros([nq 1]);
    h = zeros(nq);

    % --- Add Grad/Hess of the likelihood (matching) term
    [g1, h1] = obj.gradHessMatchingAffine();
    g = g + g1;
    h = h + h1;
    clear g1 h1
    
    % --- 2d case
    if size(obj.mu, 3) == 1
        nq = size(g, 1);
%         sub = [1 2 4 7 8 10];
        sub = [3 5 6 9 11 12];
        sub = sub(sub <= nq);
        g(sub) = 0;
        h(sub,sub) = 0;
    end

    % --- Add Grad/Hess of the prior term
    [g1, h1] = obj.gradHessPriorAffine();
    g = g + g1;
    h = h + h1;
    clear g1 h1
    
    % --- 2d case
    if size(obj.mu, 3) == 1
        diagh = diag(h);
        diagh(sub) = 1;
        h(1:nq+1:end) = diagh;
    end
    
    % --- Additional regularisation in case H is singular
    while rcond(h) < 1e-5
        r = 1e-7 * max(diag(h)) * eye(nq);
        h = h + r;
        clear r
    end
    
    
    if nargout == 0
        obj.hq.dim  = size(h);
        obj.hq(:)   = h(:);
        h           = obj.hq;
        obj.hq.utd  = true;
    end
    
end