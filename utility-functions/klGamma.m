function kl = klGamma(alpha, beta, alpha0, beta0)
% FORMAT kl = klGamma(alpha, beta, alpha0, beta0)
% alpha  - Alpha parameter of the posterior
% beta   - Beta parameter of the posterior
% alpha0 - Alpha parameter of the prior
% beta0  - Beta parameter of the prior
%
% Compute the KL-divergence from a prior to a posterior Gamma distribution:
%     KL( q || p ) = KL( n,B || n0, B0 )
% It is usually part of the lower bound in a variational framework.

    kl = (alpha - alpha0) * psi(alpha) ...
          + gammaln(alpha0) - gammaln(alpha) ...
          + alpha0 * (log(beta) - log(beta0)) ...
          + alpha * (beta0 / beta - 1);
    
end