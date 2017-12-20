function ll = lbPrecisionMatrix(EA, N, n0, EA0)
% FORMAT ll = lbPrecisionMatrix(EA, N, (n0), (EA0))
% EA  - Posterior expected precision matrix E[Az] = (n0+N)/n0 * B
% N   - Number of observations
% n0  - Prior degrees of freedom [20]
% EA0 - Prior expected precision matrix [identity]
%
% Part of the lower-bound encompassing Az terms
%     -KL(q||p) = E[ln p(Az)] - E[ln q(Az)]

    % Default parameters;
    if nargin < 4
        EA0 = eye(size(EA));
        if nargin < 3
            n0 = 20;
        end
    end

    % If n == 0 : No Wishart prior
    if n0 == 0
        ll = 0;
        return
    end
    
    % Use usual Wishart parameters
    n  = n0+N;
    B  = EA/n;
    B0 = EA0/n0;
    
    ll = -klWishart(n, B, n0, B0);
end