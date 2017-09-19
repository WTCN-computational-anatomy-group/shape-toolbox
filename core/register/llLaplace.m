function ll = llLaplace(varargin)
% FORMAT ll = llLaplace((H), (vs), (prm), ...)
%
% ** Keyword arguments **
% H   - Hessian of the negative log-likelihood w.r.t. the latent variable.
%       Its inverse is the covariance of Laplace approximation of 
%       p(Z | F, THETA).
%       If not a square matrix, it should be a (eventually symmetric)
%       tensor field.
% vs  - Voxel size of the lattice when part of the Hessian is made of L'L.
% prm - Parameters of the differential operator when part of the Hessian is
%       made of L'L.
% ** Output **
% ll - Parts of the log-likelihood of a Laplace approximation.
%
% This function returns parts of the log-likelihood of a zero-mean Gaussian 
% of precision (i.e., inverse covariance) matrix H. It is used to compute
% the log-likelihood of the Laplace approximation of a given posterior:
% log p(z | F, THETA) = -N/2 log(2pi) - 1/2 log(|H|) - 1/2 z'Hz
%
% This function returns only -N/2 log(2pi) - 1/2 log(|H|). The matching 
% part (- 1/2 z'Hz) can be computed with specfic functions 
% (llMatchingNormal, llMatchingLaplace, ...)

    % --- Parse inputs
    p = inputParser;
    p.FunctionName = 'llLaplace';
    p.addOptional('H',     [], @(X) size(X, 2) > 1 || length(X) > 2);
    p.addOptional('vs',    [], @(X) size(X, 1) == 1 && length(X) == 3);
    p.addOptional('prm',   [], @(X) size(X, 1) == 1 && length(X) == 5);
    p.addParameter('debug',  false, @isscalar);
    p.parse(varargin{:});
    H     = p.Results.H;
    vs    = p.Results.vs;
    prm   = p.Results.prm;
    debug = p.Results.debug;
    
    if debug, fprintf('* llLaplace\n'); end;
    
    ll = 0;
    
    if ~isempty(H)
        % --- Square matrix case
        if length(size(H)) == 2 && size(H, 1) == size(H, 2)
            N  = size(H, 1);
            ll = N * log(2*pi) + logDet(H);
            ll = -0.5 * ll;

        % --- Symmetric tensor case
        elseif length(size(H)) == 4
            dim = [size(H) 1 1];
            [~, N] = symIndices(dim(4));
            N  = N * prod(dim(1:3));
%             ll = N * log(2*pi) + logDetSymTensor(H);
            ll = N * log(2*pi) + sumall(log(pointwise3(H, 'd')+eps('single')));
            ll = -0.5 * ll;

        % --- Full tensor case
        elseif length(size(H)) == 5 && size(H, 4) == size(H, 5)
            dim = [size(H) 1 1];
            N   = prod(dim(1:4));
%             ll  = N * log(2*pi) + logDetTensor(H);
            ll  = N * log(2*pi) + sumall(log(pointwise3(H, 'd')));
            ll  = -0.5 * ll;

        else
            error('Unrecognized Hessian form')
        end
    end
    
    % --- L'L part
    if ~isempty(prm)
        dim = [size(H) 1 1];
        [~, ld] = spm_shoot_greens('kernel', dim(1:3), [vs prm]);
        ll = ll - 0.5 * ld(1);
        if isempty(H)
            ll = ll - 0.5 * ld(2) * log(2*pi);
        end
    end
    
end