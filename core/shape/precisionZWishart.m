function [A, n] = precisionZWishart(A0, n0, ezz, N)
% FORMAT A = precisionZWishart(A0, n0, ezz, N)
% A0  - Wishart prior for the precision matrix of Z
% n0  - Degrees of freem of the Wishart prior
% ezz - Posterior expectation of the second order statistics: sum_n E[z_n z_n']
% N   - Number of observations
%
% Update the posterior precision of the latent coordinates based on 
% its prior and some (expected) observations.
% 
% The model assumes
% - p(z|A) = N(z|A) where A is the precision of a Normal distribution
% - p(A) = W(A0, n0) where (A0,n0) are parameters of a Wishart distribution

    A = inv(A0) + ezz;
    n = n0 + N;
    
    % Stable inverse, as no ARD-style pruning is used
    [V,D] = eig(A);
    D     = loadDiag(D);
    A     = real(V * (D \ V') );
    A     = n*A; %  <- we actualy store E[X] = nA
    
end