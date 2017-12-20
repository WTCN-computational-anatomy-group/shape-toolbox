function [Q, iQ] = gnScalePG_old(ezz, n0, N, w)
% FORMAT [Q, iQ] = gnScalePG(ezz, n0, N, (w))
% ezz - Expected value of sum_n z_n * z_n'.
%       > Must have been orthogonalised before.
% n0  - Number of degrees of freedom of the Wishart prior
% N   - Number of observations
% w   - Weight on the WLW part of the z prior [1]
%
% Gauss-Newton optimisation of the scaling factor between LW and Z

    if nargin < 4
        w = 1;
    end
    
    K = size(ezz, 1);
    d = diag(ezz);
    
    % --- Initial state
    q = zeros(K, 1);
    Q = diag(exp(q));
    iQ = diag(exp(-q));
    E = obj(ezz, Q, iQ, n0, N, K, w);
%     fprintf('E = %g\n', E);
    
    % --- Gauss-Newton optimisation
    for gnit = 1:100
        
        prevE = E;
        
        dq = dir(q, d, n0, N, K, w);
        
        % --- Line search
        armijo = 1;
        for lsit=1:10
            nq = q + dq/armijo;
            nQ = diag(exp(nq));
            niQ = diag(exp(-nq));
            nE = obj(ezz, nQ, niQ, n0, N, K, w);
%             fprintf('try E = %g (%g)\n', nE, nE - E);
            if nE > E
                q  = nq;
                Q  = nQ;
                iQ = niQ;
                E  = nE;
                ok = true;
                break;
            else
                armijo = 2*armijo;
                ok = false;
            end
        end
        if ~ok
%             fprintf('Line search failure\n')
            break
        end
        if abs((E - prevE)/prevE) < 1E-5
%             fprintf('Convergence\n')
            break
        end
%         fprintf('E = %g\n', E);
        
    end

end

function e = obj(ezz, Q, iQ, n0, N, K, w)
% Compute the objective function    

    ezz = Q*ezz*Q;
    ww  = iQ*eye(K)*iQ;
    A   = precisionZWishart(n0, ezz, N);

%     e1 = - trace(ww);
%     e2 = - trace(ezz*A);
%     e3 = + w * N * proba('LogDet', ww);
%     e4 = - (N+n0)*K*trace(n0*A/(n0+N));
%     e5 = + (n0+N) * proba('LogDet', A/(n0+N));
%     e = 0.5 * (e1 + e2 + e3 + e4 + e5);

      e1 = - trace(ww);
      e2 = N * proba('LogDet', A+w*ww);
      e3 = - trace(ezz*A);
      e4 = - K*(n0+N) * trace(n0*A/(n0+N));
      e5 = - n0 * proba('LogDet', A/(n0+N));
      e  = 0.5 * (e1 + e2 + e3 + e4 + e5);
%     fprintf('%6g %6g %6g %6g %6g\n', e1, e2, e3, e4, e5);

end

function dq = dir(q, d, n0, N, K, w)
% Compute line search direction
% (Obtained with the symbolic toolbox)

%     dq = ( exp(2*q) ...
%            .* ...
%            ( n0 + d.*exp(2*q) ).^3 ...
%            .* ...
%            ( 2*exp(-2*q) ...
%              - 2*N*w ...
%              - (4*d.*exp(2*q)*(N + n0)) ...
%                ./ ...
%                (n0 + d.*exp(2*q)) ...
%              + (2*d.^2.*exp(4*q)*(N + n0)) ...
%                ./ ...
%                (n0 + d.*exp(2*q)).^2 ...
%              + (2*K*n0*d.*exp(2*q)*(N + n0)) ...
%                ./ ...
%                (n0 + d.*exp(2*q)).^2 ...
%            ) ...
%          ) ...
%          ./ ...
%          ( 4*( d.^3.*exp(6*q) ...
%                + n0^3 ...
%                + 3*n0^2*d.*exp(2*q) ...
%                + 3*n0*d.^2.*exp(4*q) ...
%                + 2*n0^3*d.*exp(4*q) ...
%                - K*n0^3*d.*exp(4*q) ...
%                + 2*N*n0^2*d.*exp(4*q) ...
%                + K*n0^2*d.^2.*exp(6*q) ...
%                - K*N*n0^2*d.*exp(4*q) ...
%                + K*N*n0*d.^2.*exp(6*q) ...
%              ) ...
%          );

    dq = -(exp(-2*q) - (d.*exp(2*q))./(n0 + d.*exp(2*q)) + (d.^2.*exp(4*q))./(n0 + d.*exp(2*q)).^2 - (d.*n0.*exp(2*q))./(n0 + d.*exp(2*q)) - (N*exp(-2*q).*(exp(4*q).*d.^2 + 2*w*exp(2*q).*d*n0 + exp(4*q).*d + w*n0^2))./((n0 + d.*exp(2*q)).*(d + n0*w*exp(-2*q) + 1)) + (K*d.*n0.*exp(2*q)*(N + n0))./(n0 + d.*exp(2*q)).^2)./(2*exp(-2*q) + (2*d.*exp(2*q))./(n0 + d.*exp(2*q)) - (6*d.^2.*exp(4*q))./(n0 + d.*exp(2*q)).^2 + (4*d.^3.*exp(6*q))./(n0 + d.*exp(2*q)).^3 + (4*N*d.*exp(2*q))./(n0 + d.*exp(2*q)) + (2*d.*n0.*exp(2*q))./(n0 + d.*exp(2*q)) - (2*d.^2.*n0.*exp(4*q))./(n0 + d.*exp(2*q)).^2 - (2*N*d.*(exp(4*q).*d.^2 + 2*w*exp(2*q).*d*n0 + exp(4*q).*d + w*n0^2))./((n0 + d.*exp(2*q)).^2.*(d + n0*w*exp(-2*q) + 1)) - (2*N*exp(-2*q).*(exp(4*q).*d.^2 + 2*w*exp(2*q).*d*n0 + exp(4*q).*d + w*n0^2))./((n0 + d.*exp(2*q)).*(d + n0*w*exp(-2*q) + 1)) + (4*K*d.^2.*n0.*exp(4*q)*(N + n0))./(n0 + d.*exp(2*q)).^3 - (2*K*d.*n0.*exp(2*q).*(N + n0))./(n0 + d.*exp(2*q)).^2 + (2*N*n0*w*exp(-4*q).*(exp(4*q).*d.^2 + 2*w*exp(2*q).*d*n0 + exp(4*q).*d + w*n0^2))./((n0 + d.*exp(2*q)).*(d + n0*w*exp(-2*q) + 1).^2));
    dq = -dq;
end

%% ________________________________________________________________________
%
% Symbolic code to obtain the search direction
%  --------------------------------------------
%
% Wishart Az with Identiy prior and n0 degrees of freedom
% We assume that orthogonalisation has been performed before.
%
% NB: w is the weight on the WLW part of the z prior
%
% Assumptions:
% W'L'LW = I     (Identity)
% E[Z]E[Z]' = D  (Diagonal)
% Sz ~ 0         (Neglectable)
% => E[ZZ'] ~ D  (Almost diagonal)
% B0 = I/n0      (Expected value = identity)
%
% q: K rescaling factors such that
% Q  = diag(exp(q))
% iQ = diag(exp(-q))
% => Q  = Q'
%    iQ = iQ'
%
% Diagonal simplification
% E[Az] = (n0+N)B = (n0+N)inv(iB0 + E[ZZ']) = (n0+N)inv(n0I+D)
% => E[Az]k = (n0+N)/(n0+dk)
% QE[ZZ']Q (k) = dk exp(2 qk)
% E[Az] ~ Q (k) = (n0+N)/(n0+dk exp(2 qk))
%
% Terms
% 1) -Tr(iQW'L'LWiQ) = -Tr(iQiQ) = -Tr(iQ^2) = -sum_k  exp(-2 qk)
% 2) wN logdet(iQW'L'LWiQ) = wN logdet(iQ^2) = 2wN logdet(iQ) = -2wN sum_k qk
% 3) -Tr(QE[ZZ']Q E[Az~Q]) = -sum_k (n0+N) dk exp(2 qk) / (n0+dk exp(2 qk))
% 4) -(n0+N)K Tr(iB0B) = -(n0+N)*K*n0/(n0+dk exp(2 qk))
% 5) (n0+N) logdet(B) = -(n0+N) sum_k log(n0 + dk exp(2 qk))
% _________________________________________________________________________
% syms n0 N K d q w
% E1 = -exp(-2*q);
% E2 = -w*2*N*q;
% E3 = -(n0+N)*d*exp(2*q) / (n0 + d*exp(2*q));
% E4 = -(n0+N)*n0*K/(n0 + d*exp(2*q));
% E5 = -(n0+N) * log(n0 + d * exp(2*q));
% _________________________________________________________________________
% E = simplify(E1+E2+E3+E4+E5);
% dE = simplify(diff(E, q));
% d2E = simplify(diff(dE, q));
% dir = simplify(-d2E\dE);