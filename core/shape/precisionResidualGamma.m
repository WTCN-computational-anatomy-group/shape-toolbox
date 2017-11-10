function lam = precisionResidualGamma(lam0, n0, err, N, lat)
% FORMAT lam = precisionResidualGamma(lam0, n0, err, N, lat)
% lam0  - Expected value aof the Gamma prior
% n0    - "Degrees of freedom" of the Gamma prior 
%         (alpha0 = 0.5*prod(lat)*3*n0)
% err   - sum of all individual contributions, i.e. 0.5*trace(E[RR]LL)
% N     - Number of observations (subjects)
% lat   - Dimension of the residual lattice
% lam   - Expected value of the Gamma posterior
%
% Update the posterior precision of the residual field based on 
% its prior and some (expected) observations.

    alpha0 = 0.5*prod(lat)*3*n0;
    beta0  = alpha0/lam0;
    alpha  = alpha0 + 0.5*prod(lat)*3*N;
    beta   = beta0  + err;
    lam    = alpha/beta;
    lam    = double(lam);
    
end