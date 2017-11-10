function kl = klWishart(n, B, n0, B0)
% FORMAT kl = klWishart(n, B, n0, B0)
% n  - DF parameter of the posterior
% B  - Precision parameter of the posterior
% n0 - DF parameter of the prior
% B0 - Precision parameter of the prior
%
% Compute the KL-divergence from a prior to a posterior Wishart 
% distribution:
%     KL( q || p ) = KL( n,B || n0, B0 )
% It is usually part of the lower bound in a variational framework.
    
    p = size(B, 1);
    
    kl =   0.5 * n * (trace(B0\B) - p) ...
         + proba('LogGamma', n0/2, p) - proba('LogGamma', n/2, p) ...
         + 0.5 * (n - n0) * proba('DiGamma', n/2, p) ...
         + 0.5 * n0 * ( proba('LogDet', B0) - proba('LogDet', B) );

end