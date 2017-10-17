function A = precisionZWishart(n0, ezz, N, md)
% FORMAT A = precisionZWishart(n0, ezz, N)
% n0  - Degrees of freedom of the Wishart prior
% ezz - Posterior expectation of the second order statistics: sum_n E[z_n z_n']
% N   - Number of observations
% A   - Expected value of the Wishart posterior
% md  - 'mean' [default] or 'mode'
%
% Update the posterior precision of the latent coordinates based on 
% its prior and some (expected) observations.
% 
% The model assumes
% - p(z|A) = N(z|A) where A is the precision of a Normal distribution
% - p(A) = W(A0, n0) where (A0,n0) are parameters of a Wishart distribution

    if nargin < 4
        md = 'mean';
    end

    A = n0*eye(size(ezz)) + ezz;
    n = n0 + N;
    
    % Stable inverse, as no ARD-style pruning is used
    [V,D] = eig(A);
    D     = loadDiag(D);
    A     = real(V * (D \ V') );
    
    if strcmpi(md, 'mean')
        A = n*A; %  <- we actualy store E[X] = nA
    elseif strcmpi(md, 'mode')
        A = (n + size(ezz, 1) - 1)*A;
    else
        error('Unknown mode %s', md);
    end
    
end