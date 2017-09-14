function ll = llLaplace(varargin)
% FORMAT ll = llLaplace(('latent', hz), ('residual', hr), ('affine', hq))
%
% ** Keyword arguments **
% hz - Hessian of the log-likelihood w.r.t. current Z.
%        Its inverse is the covariance of Laplace approximation of 
%        p(Z | F, Mu, W, R).
% hr - Hessian of the log-likelihood w.r.t. current R.
%        Its inverse is the covariance of Laplace approximation of 
%        p(R | F, Mu, W, Z).
% hq - Hessian of the log-likelihood w.r.t. current Q.
%        Its inverse is the covariance of Laplace approximation of 
%        p(q | F, Mu, W, Z).
% ** Output **
% ll - Parts of the log-likelihood which encompass Laplace approximation
%        terms.
%
% Parts of Laplace approximation of the posterior:
% p(F|Mu) = int_{Z,R} p(F | Mu, Z, R) dZ dR
%         = (2pi)^(-K*N*Q/2) |-H(R*,Z*,q*)|^(1/2) p(F | Mu, Z*, R*, q*) p(Z*) p(R*) p(q*)
%         = (2pi)^(-K*N*Q/2) |-H(R*)|^(1/2) |-H(Z*)|^(1/2) |-H(Q*)|^(1/2) p(F | Mu, Z*, R*, q*) p(Z*) p(R*) p(q*)
% where Z*, R*, q* define a mode of p(F|Mu,Z,R,Q) and H is the hessian of
% respectively {R,Z,Q}, Z, R qnd Q around this mode.
% This function thus returns:
% ll = -0.5 * KNQ * log(2pi) + 0.5 * log(|-H(Z*)|) + 0.5 * log(|-H(R*)|) + 0.5 * log(|-H(q*)|)

    % --- Parse inputs
    p = inputParser;
    p.FunctionName = 'llLaplace';
    p.addParameter('latent',    [], @checkarray);
    p.addParameter('residual',  [], @checkarray);
    p.addParameter('affine',    [], @checkarray);
    p.addParameter('debug',  false, @isscalar);
    p.parse(varargin{:});
    hz = p.Results.latent;
    hr = p.Results.residual;
    hq = p.Results.affine;
    
    if p.Results.debug, fprintf('* llLaplace\n'); end;
    
    ll = 0;
    
    % --- Z part
    if ~isempty(hz)
        K  = size(hz, 1);
        ll = ll - logDet(hz);
        useZ = 1;
    else
        K = 1;
        useZ = 0;
    end
    
    % --- R part
    if ~isempty(hr)
        dim = [size(hr) 1 1 1];
        N   = prod(dim(1:4));
        ll  = ll - logDetSymTensor(hr);
        useR = 1;
    else
        N = 1;
        useR = 0;
    end
    
    % --- Q part
    if ~isempty(hq)
        dim = [size(hq) 1 1 1];
        Q   = prod(dim(1:4));
        ll  = ll - logDet(hq);
        useQ = 1;
    else
        Q = 1;
        useQ = 0;
    end
    
    % --- Common part
    ll = ll * (-1)^(N*useR + K*useZ + Q*useQ);
    ll = ll + N^useR * K^useZ * Q^useQ * log(2 * pi);
    
    % --- Common factor
    ll = -0.5 * ll;
    
end