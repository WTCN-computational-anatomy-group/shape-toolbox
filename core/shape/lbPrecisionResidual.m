function ll = lbPrecisionResidual(Elam, N, n0, Elam0, lat)
% FORMAT ll = lbPrecisionResidual(Elam, N, n, Elam0, lat)
% Elam  - Posterior expected precision E[lam] = (n0+N)/n0 * lam0
% N     - Number of observations
% n0    - Prior degrees of freedom
% Elam0 - Prior expected precision matrix
% lat   - Dimension of the velocity lattice
%
% Part of the lower-bound encompassing Az terms
%     -KL(q||p) = E[ln p(Az)] - E[ln q(Az)]

    % Default parameters;
    if nargin < 4
        Elam0 = 10;
        if nargin < 3
            n0 = 20;
        end
    end

    % Use usual Wishart parameters
    alpha0 = 0.5 * prod(lat) * 3 * n0;
    alpha  = alpha0 + 0.5 * prod(lat) * 3 * N;
    beta   = alpha/Elam;
    beta0  = alpha0/Elam0;
    
    ll = -spm_prob('Gamma', 'kl', alpha, beta, alpha0, beta0);
end