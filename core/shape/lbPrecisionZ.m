function ll = lbPrecisionZ(EA, N, n0)
% FORMAT ll = lbPrecisionZ(EA, N, (n0), (B0))
% EA - Posterior expected precision matrix E[Az] = (n0+N)B
% N  - Number of observations
% n0 - Prior degrees of freedom [20]
% B0 - Prior expected precision matrix [identity]
%
% Part of the lower-bound encompassing Az terms
% > E[ln p(Az)] - E[ln q(Az)]

    % Default parameters;
    if nargin < 3
        n0 = 20;
    end

    % If n == 0 : No Wishart prior
    if n0 == 0
        ll = 0;
        return
    end
    
    % Use Wikipedia's definition of the parameters
    B  = EA/(n0+N);
    K  = size(EA, 1);

    % Compute lower bound
    ll = - N/2 * proba('DiGamma', (N+n0)/2, K) ...
        + proba('LogGamma', (N+n0)/2, K) ...
        - proba('LogGamma', n0/2, K) ...
        + (n0+N)*K/2 * (1 - trace(n0*B)) ...
        + n0/2*( proba('LogDet', B) - proba('LogDet', eye(K)/n0));
end