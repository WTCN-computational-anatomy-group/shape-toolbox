function kl = klNormal(mu, A, mu0, A0, mode)
% FORMAT kl = klNormal(mu, A, mu0, A0, mode)
% mu   - Mean parameter of the posterior
% A    - Covariance/precision parameter of the posterior
% mu0  - Mean parameter of the prior
% A0   - Covariance/precision parameter of the prior
% mode - 'covariance' or 'precision'
%
% Compute the KL-divergence from a prior to a posterior Normal 
% distribution:
%     KL( q || p ) = KL( mu, A || mu0, A0 )
% It is usually part of the lower bound in a variational framework.
    
    if nargin < 5
        mode = 'covariance';
    end

    k = size(A, 1);
    x = mu0-mu;
    if strcmpi(mode, 'covariance')
        kl = trace(A0\A) + x(:)' * (A0 \ x(:)) ...
             + proba('LogDet', A0) - proba('LogDet', A) - k;
    else
        kl = trace(A\A0) + x(:)' * A0 * x(:) ...
             + proba('LogDet', A) - proba('LogDet', A0) - k;
    end
    
    kl = 2 * kl;

end