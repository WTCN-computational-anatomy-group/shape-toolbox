function ld = ldapprox(varargin)
% _________________________________________________________________________
%
%      Stochastic approximation of logDet(H+L)
%
% -------------------------------------------------------------------------
% FORMAT ld = ldapprox(L, ...)
% FORMAT ld = ldapprox(prm, ...)
% FORMAT ld = ldapprox(prm, H, ...)
%
% MANDATORY
% ---------
% L       - Spare or full matrix
%    or
% prm     - Parameters of the prior precision (if spm_diffeo/spm_field)
%
% OPTIONAL
% --------
% H       - Hessian of the data term:         nx*ny*nz*6 or [empty/0]
%
% KEYWORD ARGUMENTS
% -----------------
% vs      - Lattice voxel size                              [1 1 1]
% fmg     - Full Multigrid options                          [cyc rlx]
%           cyc - Number of Full Multigrid cycles           [2]
%           rlx - Number of relaxation iterations per cycle [2]
% moments - Number of moments                               [5]
% samples - Number of random samples for moment estimation  [5]
% mc      - Number of samples for Monte Carlo integration   [10]
% method  - Sampling method                      'gaussian'/['rademacher']
% type    - Input type                              'field'/['diffeo']
%
% OUTPUT
% ------
% ld       - Log-determinant
%
% -------------------------------------------------------------------------
% Granziol, Roberts & Osborne
% "VBALD - Variational Bayesian Approximation of Log Determinants"
% https://arxiv.org/abs/1802.08054
% _________________________________________________________________________

    % Note: in their example on sparse matrices, Granziol et al. use 5
    % moments and 5 samples. They might be reasonable default values.

    % ---------------------------------------------------------------------
    % Parse inputs
    % ---------------------------------------------------------------------
    p = inputParser;
    p.addRequired('L');
    p.addOptional('H',        []);
    p.addParameter('vs',      [1 1 1]);
    p.addParameter('fmg',     [2 2]);
    p.addParameter('moments', 5);
    p.addParameter('samples', 5);
    p.addParameter('mc',      10);
    p.addParameter('method',  'rademacher');
    p.addParameter('type',    'diffeo');
    p.addParameter('dim',     []);
    p.parse(varargin{:});
    
    L     = p.Results.L;
    H     = p.Results.H;
    vs    = p.Results.vs;
    fmg   = p.Results.fmg;
    nm    = p.Results.moments;
    ns    = p.Results.samples;
    mc    = p.Results.mc;
    samp  = p.Results.method;
    type  = p.Results.type;
    dim  = p.Results.dim;

    % matrix or paramter form?
    if size(L,1) ~= size(L,2)
        prm = L;
    end
    
    % empty Hessian?
    if isempty(H) || numel(H) == 1
        if isempty(dim)
            lat = [size(L, 1) 1];
            nc  = 1;
        end
        ulam = 0;
    else
        dimh = [size(H) 1 1];
        lat = dimh(1:3);
        [~,nc]  = spm_matcomp('SymIndices', dimh(4), 'k');
        H = single(H);
        ulam = max(max(max(sum(abs(H),4))));
    end
    
    % Compute upper bound on eigenvalues
    if ~isempty(L)
        ulam = ulam + sum(abs(full(L(1,:))));
    else
        ulam = ulam + sumall(abs(spm_diffeo('kernel',lat,[vs prm])));
    end
    
    % Compute moments
    mom = trapprox(L/ulam, H/ulam, 'vs', vs, 'fmg', fmg, ...
                   'moments', nm, 'samples', ns, 'matrix', 'H+L', ...
                   'method', samp, 'type', type, 'dim', dim);
    mom = mom/prod([lat nc]);
                
    % Compute beta parameters (Maximum Likelihood)
    alpha = mom(1)*(mom(1)-mom(2))/(mom(2)-mom(1)^2);
    beta  = (1/mom(1)-1)*alpha;
    
    if alpha > 0 && beta > 0    % We can use the beta pior
        prior = 'beta';
    else                        % We have to use the uniform prior
        prior = 'uni';
    end
    
    % Compute coefficients
    coef = VBALD(mom, 1e-3, mc, prior, alpha, beta);

    % logdet(H+L) = K * ( E[log(lam)] + log(ulam) )
    K = prod([lat nc]);
    ld = K * ( mcmc_logq(coef, mc, prior, alpha, beta) + log(ulam) );
        
end

function [c,l] = VBALD(m,tol,ns,prior,a,b)
% FORMAT c = vbald(m, tol)
% m    - moments [K 1]
% tol  - tolerance
% c    - coefficients [K 1]
% l    - objective function
%
% Compute power coefficients by Gauss-Newton optimisation of 

    c0   = zeros(size(m));
    gtol = inf;
    l0   = mcmc_gh(c0,ns,prior,a,b) + sum(c0.*m);
    i    = 1;
%     fprintf('i = %3i ;         ;            ; l = %6f\n', 0, l0);
    armijo = 1;
    while gtol > tol
        [~,g,h] = mcmc_gh(c0,ns,prior,a,b);
        g = m - g;
        d = -h\g;
        ok     = false;
        for it=1:20
            c = c0 - d/armijo;
            l = mcmc_gh(c,ns,prior,a,b) + sum(c.*m);
%             fprintf('i = %3i ; j = %3i ; a = %9.3g ; l = %6f\n', i, it, armijo, l);
            if l < l0
                armijo = 0.9 * armijo;
                ok = true;
                break
            end
            armijo = armijo*2;
        end
        if ok
            gtol = abs(l0-l);
            l0 = l;
            c0 = c;
            i  = i+1;
        else
            l = l0;
            c = c0;
            return
        end
    end
end

function l = q(lam,c,prior,a,b)
    switch prior(1)
        case 'u'
            l = q_uni(lam,c);
        case 'b'
            l = q_beta(lam,c,a,b);
    end
end

function l = q_beta(lam,c,a,b)
    l = betapdf(lam,a,b) * factexp(lam,c);
end

function l = q_uni(lam,c)
    l = factexp(lam,c);
end

function e = factexp(lam,c)
    e = exp(-( 1 + sum(c(:).*(lam.^(1:length(c))')) ));
end

function [l,g,h] = mcmc_gh(c,ns,prior,a,b)
% Compute E[q(lam)*lam^j] by Monte Carlo integration
% a,b - Parameters of the prior Beta distribution
% c   - Power coefficients
% j   - Power in the integral
% ns  - Number of samples
    switch nargout
        case 1
            np = 1;
        case 2
            np = 1 + length(c);
        case 3
            np = 1 + 2*length(c);
    end
    
    s = zeros(np,1);
    for i=1:ns
        switch prior(1)
            case 'u'
                lam = unifrnd(eps,1);
            case 'b'
                lam = betarnd(a,b);
        end
        cum = q(lam,c,prior,a,b);
        s(1) = cum;
        cum = log(cum);
        for j=2:length(s)
            cum = cum + log(lam);
            s(j) = s(j) + exp(cum);
        end
    end
    s = s/ns;
    
    l = s(1);
    if nargout > 1
        g = s(2:length(c)+1);
        g = g(:);
        if nargout > 2
            h = zeros(length(c));
            for j=1:length(c)
                for k=j+1:length(c)
                    h(j,k) = s(1+j+k);
                    h(k,j) = s(1+j+k);
                end
                h(j,j)= s(1+j+j);
            end
%             h = spm_matcomp('LoadDiag', h);
            h = h + eye(size(h)) * max(diag(h)) * 1e-7;
        end
    end
end


function s = mcmc_logq(c,ns,prior,a,b)
% Compute E[log(lam)*exp(-(1+sum c_j*lam^j))] by Monte Carlo integration
% a,b - Parameters of the prior Beta distribution
% c   - Power coefficients
% ns  - Number of samples
    s = 0;
    for i=1:ns
        switch prior(1)
            case 'u'
                lam = unifrnd(eps,1);
            case 'b'
                lam = betarnd(a,b);
        end
        s = s + log(lam) * factexp(lam,c);
    end
    s = s/ns;
end