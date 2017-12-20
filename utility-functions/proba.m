function varargout = proba(id, varargin)
% Collection of tools for dealing with probability distributions.
%
% FORMAT ld = proba('LogDet', A)
% > Robust log of determinant of positive-definite matrix
%
% FORMAT lg = proba('LogGamma', a, p)
% > Log of multivariate gamma function of order p
%
% FORMAT dg = proba('DiGamma', a, p)
% > Multivariate digamma function of order p
%
% FORMAT ll = proba('ELogDetWishart', B, n)
% > Expected log-determinant of a Wishart distribution.
%
% FORMAT ll = proba('ELogDetNonCentralWishart', MM, n, bound)
% > Lower or upper bound of the expected log-determinant of a 
%   non-central Wishart distribution.
%
% FORMAT ld = proba('LogDetDiffeo', lat, vs, prm)
% > Log determinant of the prior precision on velocities

    switch lower(id)
        case 'logdet'
            [varargout{1:nargout}] = LogDet(varargin{:});
        case 'loggamma'
            [varargout{1:nargout}] = LogGamma(varargin{:});
        case 'digamma'
            [varargout{1:nargout}] = DiGamma(varargin{:});
        case 'elogdetwishart'
            [varargout{1:nargout}] = ELogDetWishart(varargin{:});
        case 'elogdetnoncentralwishart'
            [varargout{1:nargout}] = ELogDetNonCentralWishart(varargin{:});
        case 'logdetdiffeo'
            [varargout{1:nargout}] = LogDetDiffeo(varargin{:});
        otherwise
            error('Unknown function %s', id);
    end
end

% -------------------------------------------------------------------------
function ld = LogDet(A)
% FORMAT ld = LogDet(A)
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

% -------------------------------------------------------------------------
function lg = LogGamma(a, p)
    if nargin < 2
        p = 1;
    end
    lg = (p*(p-1)/4)*log(pi);
    for i=1:p
        lg = lg + gammaln(a + (1-p)/2);
    end
end

% -------------------------------------------------------------------------
function dg = DiGamma(a, p)
    if nargin < 2
        p = 1;
    end
    dg = 0;
    for i=1:p
        dg = dg + psi(a + (1-p)/2);
    end
end

% -------------------------------------------------------------------------
function ll = ELogDetWishart(B, n)
% FORMAT ll = proba('ELogDetWishart', B, n)
%
% Expectation of the log determinant of a precision matrix that stems from
% a Wishart distribution.

    K = size(B,1);
    % K ln 2 + log|B| + psi_K(n/2)
    ll = K * log(2) + LogDet(B) + DiGamma(n/2, K);
end

% -------------------------------------------------------------------------
function ll = ELogDetNonCentralWishart(MM, n, bnd)
% FORMAT ll = proba('ELogDetNonCentralWishart', MM, n, bnd)
%
% MM  - Non centrality matrix M'M
% n   - degrees of freedom
% bnd - 'lower' [default] or 'upper'
%
% Lower or upper bound of the Expected value of log-determinant of a 
% non-central Wishart distribution of parameters (B, MM, n).
% Y stems from a non-central Wishart distribution if
% Y = (X+M)'(X+M), with
%   - X is p*n and Xi ~ N(0,I_p)
%   - M is fixed p*n
%   - n > p - 1
%
% Upper and lower bounds were dervived by Kim and Lapidoth in:
% "On the log determinant of noncentral Wishart matrices", 
% Proceedings of the IEEE International Symposium on Information Theory,
% 2003. 

    if nargin < 4
        bnd = 'lower';
    end
    
    S = svd(MM);
    p = size(MM);
    
    switch lower(bnd)
        case 'upper'
            ll = 0;
            for k=1:p
                ll = ll + gn(S(k), n);
            end
        
        case 'lower'
            ll = gn(S(1), n);
            for k=2:p
                ll = ll + psi(n-k+1);
            end
            
        otherwise
            error('Unknwon bound %s', bnd)
    end

end

function ll = gn(lambda, n)
    ll = log(lambda) - expint(-lambda);
    for i=1:(n-1)
        ll = ll + (-1/lambda)^i * ...
                  ( exp(-lambda)*factorial(i-1) - ...
                    factorial(n-1)/(i*factorial(n-1-i)) );
    end
end

% -------------------------------------------------------------------------
function ld = LogDetDiffeo(lat, vs, prm)
    [~, ld] = spm_shoot_greens('kernel', lat, [vs prm]);
    ld = ld(1);
end