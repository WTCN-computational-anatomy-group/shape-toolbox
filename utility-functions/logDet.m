function ld = logDet(A)
% FORMAT ld = logDet(A)
% A  - A square matrix
% ld - Logarithm of determinant of A
%
% Log-determinant of a matrix
% Cholesky factorisation is used to compute a more stable log-determinant.
%__________________________________________________________________________
% Copyright (C) 2017 Wellcome Trust Centre for Neuroimaging

% John Ashburner
% $Id$

    % Cholseki decomposition of A (A = C' * C, with C upper-triangular)
    [C, p] = chol(numeric(A));
    
    if p > 0
       % A should usually be positive definite, but check anyway.
       warning(['Attempting to compute log determinant of matrix ' ...
                'that is not positive definite (p=%d).'], p);
    end

    % Because C is triangular, |C| = prod(diag(C))
    % Hence: log|C| = sum(log(diag(C)))
    % And:   log|A| = log|C'*C| = log(|C|^2) = 2 * sum(log(diag(C)))
    ld = 2 * sum(log(diag(C)));

end