function kl = klGamma(alpha, beta, alpha0, beta0, mode, K)
% _________________________________________________________________________
% FORMAT kl = klGamma(alpha, beta, alpha0, beta0)
% alpha  - Posterior alpha parameter
% beta   - Posterior beta parameter
% alpha0 - Prior alpha parameter
% beta0  - Prior beta parameter
% _________________________________________________________________________
% FORMAT kl = klGamma(n, lam, n0, lam0, 'precision', (K))
% n    - Posterior degrees of freedom
% lam  - Posterior precision
% n0   - Prior degrees of freedom
% lam0 - Prior precision
% K    - Space dimension (Multivariate case) [1]
% _________________________________________________________________________
% Compute the KL-divergence from a prior to a posterior Gamma distribution:
%     KL( q || p ) = KL( a,b || a0, b0 )
% It is usually part of the lower bound in a variational framework.
% _________________________________________________________________________

    if nargin < 6
        K = 1;
        if nargin < 5
            mode = 'common';
        end
    end
    precision = strcmpi(mode, 'precision');

    if precision
        kl = klGamma(K*alpha/2,  K*alpha/(2*beta), ...
                     K*alpha0/2, K*alpha0/(2*beta0));
    else
        kl = (alpha - alpha0) * psi(alpha) ...
              + gammaln(alpha0) - gammaln(alpha) ...
              + alpha0 * (log(beta) - log(beta0)) ...
              + alpha * (beta0 / beta - 1);
    end
    
end