function [Q, iQ] = scalePG(N, K)
% FORMAT [Q, iQ] = scalePG(N, K)
%
% Precision matrices scaling when there is no Wishart prior on Az

    % NB: This is if we use p(z|W,A) = p(z|W)p(z|A), without true
    % normalisation (i.e., p(z|W,A) is not really a density)
    q = (-log(2*N)/2) * eye(K);
    
    % If we use a normalisation, i.e., p(z|W,A) \propto p(z|W)p(z|A),
    % we should use
    % > q = (-log(N)/2) * eye(K);

    Q = diag(q);
    iQ = diag(1./q);

end