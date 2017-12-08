function kl = klWishart(n, B, n0, B0, mode)
% _________________________________________________________________________
% FORMAT kl = klWishart(n, B, n0, B0)
% n    - Posterior DF parameter
% B    - Posterior scale matrix parameter
% n0   - Prior DF parameter
% B0   - Prior scale matrix parameter
% _________________________________________________________________________
% FORMAT kl = klWishart(n, A, n0, A0, 'expectation')
% n    - Posterior DF parameter
% A    - Posterior precision parameter
% n0   - Prior DF parameter
% A0   - Prior precision parameter
% _________________________________________________________________________
% Compute the KL-divergence from a prior to a posterior Wishart 
% distribution:
%     KL( q || p ) = KL( n,B || n0, B0 )
% It is usually part of the lower bound in a variational framework.
% _________________________________________________________________________
    
    if nargin < 5
        mode = 'scale';
    end
    
    if strcmpi(mode, 'expectation')
        kl = klWishart(n, B/n, n0, B/n0);
    else
        p = size(B, 1);

        kl =   0.5 * n * (trace(B0\B) - p) ...
             + proba('LogGamma', n0/2, p) - proba('LogGamma', n/2, p) ...
             + 0.5 * (n - n0) * proba('DiGamma', n/2, p) ...
             + 0.5 * n0 * ( proba('LogDet', B0) - proba('LogDet', B) );
    end
end