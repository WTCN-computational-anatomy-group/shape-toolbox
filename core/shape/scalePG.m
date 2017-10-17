function [Q, iQ] = scalePG(N, K)
% FORMAT [Q, iQ] = scalePG(N, K)
%
% Precision matrices scaling when there is no Wishart prior on Az

    q = -log(N)/2 * ones(K,1);

    Q = diag(exp(q));
    iQ = diag(exp(-q));

end