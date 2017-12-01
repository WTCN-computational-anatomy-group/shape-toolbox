function [Q, iQ] = gnScalePG_vel(Sz, zz, ww, n0, N)

    K  = size(Sz, 1);
    s  = diag(Sz);
    z  = diag(zz);
    w  = diag(ww);
    a0 = ones(K,1);
    
    % --- Initial state
    q = zeros(K, 1);
    Q = diag(exp(q));
    iQ = diag(exp(-q));
    E = objective(Q, iQ, Sz, zz, ww, eye(K), n0, N);
%     fprintf('E = %g\n', E);
    
    % --- Gauss-Newton optimisation
    for gnit = 1:100
        
        prevE = E;
        
        dq = dir(q, s, z, w, a0, n0, N);
        if max(dq) > 7
            dq = dq / max(dq) * 7;
        end
        if min(dq) < -7
            dq = dq / min(dq) * (-7);
        end
        
        % --- Line search
        armijo = 1;
        for lsit=1:10
            nq = q + dq/armijo;
            nQ = diag(exp(nq));
            niQ = diag(exp(-nq));
            nE = objective(nQ, niQ, Sz, zz, ww, eye(K), n0, N);
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

function dq = dir(q, s, z, w, a0, n0, N, pos)
% Compute line search direction
% (Obtained with the symbolic toolbox)

%    dq = -(2.*s.*exp(2.*q) - 2.*w.*exp(-2.*q) - (2.*a0.*exp(2.*q).*(N + n0).^2.*(s + z))./(n0 + a0.*exp(2.*q).*(s + z)).^2 - (2.*a0.^2.*exp(4.*q).*(N + n0).*(s + z).^2)./(n0 + a0.*exp(2.*q).*(s + z)).^2 + (2.*a0.*exp(2.*q).*(N + n0).*(s + z))./(n0 + a0.*exp(2.*q).*(s + z)) + (2.*a0.^2.*exp(2.*q).*(N + n0).*(s + z).*(n0 + 1))./(n0 + a0.*exp(2.*q).*(s + z)).^2)./(4.*s.*exp(2.*q) + 4.*w.*exp(-2.*q) - (4.*a0.*exp(2.*q).*(N + n0).^2.*(s + z))./(n0 + a0.*exp(2.*q).*(s + z)).^2 - (12.*a0.^2.*exp(4.*q).*(N + n0).*(s + z).^2)./(n0 + a0.*exp(2.*q).*(s + z)).^2 + (8.*a0.^3.*exp(6.*q).*(N + n0).*(s + z).^3)./(n0 + a0.*exp(2.*q).*(s + z)).^3 + (4.*a0.*exp(2.*q).*(N + n0).*(s + z))./(n0 + a0.*exp(2.*q).*(s + z)) + (8.*a0.^2.*exp(4.*q).*(N + n0).^2.*(s + z).^2)./(n0 + a0.*exp(2.*q).*(s + z)).^3 - (8.*a0.^3.*exp(4.*q).*(N + n0).*(s + z).^2.*(n0 + 1))./(n0 + a0.*exp(2.*q).*(s + z)).^3 + (4.*a0.^2.*exp(2.*q).*(N + n0).*(s + z).*(n0 + 1))./(n0 + a0.*exp(2.*q).*(s + z)).^2);
    dq = -(2.*s.*exp(2.*q) - 2.*w.*exp(-2.*q) + (2.*exp(2.*q).*(N + n0).*(s + z))./(n0./a0 + exp(2.*q).*(s + z)) - (2.*exp(4.*q).*(N + n0).*(s + z).^2)./(n0./a0 + exp(2.*q).*(s + z)).^2 - (2.*a0.*exp(2.*q).*(N + n0).^2.*(s + z))./(n0 + a0.*exp(2.*q).*(s + z)).^2 + (2.*exp(2.*q).*(N + n0).*(s + z).*(n0 + 1))./(n0./a0 + exp(2.*q).*(s + z)).^2)./(((4.*s.*exp(2.*q) + 4.*w.*exp(-2.*q) + (4.*exp(2.*q).*(N + n0).*(s + z))./(n0./a0 + exp(2.*q).*(s + z)) - (12.*exp(4.*q).*(N + n0).*(s + z).^2)./(n0./a0 + exp(2.*q).*(s + z)).^2 + (8.*exp(6.*q).*(N + n0).*(s + z).^3)./(n0./a0 + exp(2.*q).*(s + z)).^3 - (4.*a0.*exp(2.*q).*(N + n0).^2.*(s + z))./(n0 + a0.*exp(2.*q).*(s + z)).^2 - (8.*exp(4.*q).*(N + n0).*(s + z).^2.*(n0 + 1))./(n0./a0 + exp(2.*q).*(s + z)).^3 + (8.*a0.^2.*exp(4.*q).*(N + n0).^2.*(s + z).^2)./(n0 + a0.*exp(2.*q).*(s + z)).^3 + (4.*exp(2.*q).*(N + n0).*(s + z).*(n0 + 1))./(n0./a0 + exp(2.*q).*(s + z)).^2)./(s.*exp(2.*q) + w.*exp(-2.*q) + (N + n0).^2./(n0 + a0.*exp(2.*q).*(s + z)) - ((N + n0).*(n0 + 1))./(n0./a0 + exp(2.*q).*(s + z)) + (exp(2.*q).*(N + n0).*(s + z))./(n0./a0 + exp(2.*q).*(s + z))) - (2.*s.*exp(2.*q) - 2.*w.*exp(-2.*q) + (2.*exp(2.*q).*(N + n0).*(s + z))./(n0./a0 + exp(2.*q).*(s + z)) - (2.*exp(4.*q).*(N + n0).*(s + z).^2)./(n0./a0 + exp(2.*q).*(s + z)).^2 - (2.*a0.*exp(2.*q).*(N + n0).^2.*(s + z))./(n0 + a0.*exp(2.*q).*(s + z)).^2 + (2.*exp(2.*q).*(N + n0).*(s + z).*(n0 + 1))./(n0./a0 + exp(2.*q).*(s + z)).^2).^2./(s.*exp(2.*q) + w.*exp(-2.*q) + (N + n0).^2./(n0 + a0.*exp(2.*q).*(s + z)) - ((N + n0).*(n0 + 1))./(n0./a0 + exp(2.*q).*(s + z)) + (exp(2.*q).*(N + n0).*(s + z))./(n0./a0 + exp(2.*q).*(s + z))).^2).*(s.*exp(2.*q) + w.*exp(-2.*q) + (N + n0).^2./(n0 + a0.*exp(2.*q).*(s + z)) - ((N + n0).*(n0 + 1))./(n0./a0 + exp(2.*q).*(s + z)) + (exp(2.*q).*(N + n0).*(s + z))./(n0./a0 + exp(2.*q).*(s + z))));

end

% function lb = objective(q, s, z, w, a0, n0, N)
function lb = objective(Q, iQ, Sz, zz, WW, A0, n0, N)
% Compute (restricted) lower bound
% (Obtained with the symbolic toolbox)

%    lb = s - w.*exp(-2.*q) + (a0.*(N + n0).*(n0 + 1))./(n0 + a0.*exp(2.*q).*(s + z)) + (exp(2.*q).*(N + n0).^2.*(s + z))./(n0 + a0.*exp(2.*q).*(s + z)) - (a0.*exp(2.*q).*(N + n0).*(s + z))./(n0 + a0.*exp(2.*q).*(s + z));

    Sz = Q*Sz*Q';
    zz = Q*zz*Q';
    WW = iQ'*WW*iQ;
    A = (n0 + N) * inv(n0 * inv(A0) + Sz + zz);

    lb = - trace((Sz+zz)*A)  ...
         + log(det(Sz)) ...
         + (n0+1)*proba('LogDet',A) ...
         + (n0+N)*trace(inv(A0)*A) ...
         - trace(WW);
     lb = 0.5 * lb;

end

%%

% -------------------------------------------------------------------------
% Symbolic solving
% -------------------------------------------------------------------------
%
% Let us write:
%   S         : Sum of individual covariance matrices of Z (= cov[Z])
%   C         : E[ZZ'] = cov[Z] + E[Z]E[Z]'
%   WW        : Covariance steming from the PGs
%   A, (N+n0) : Posterior covariance of Z
%   A0, n0    : Prior covariance of Z
%
% We want to maximise w.r.t. to a transform:
%   2*L = - trace(C*A)                      | -KL(z)
%         - log(det(S))                     | -KL(z)
%         + (n0+1)log(det(A))               | -KL(A) (n0) - KL(z) (1)
%         - (n0+N)*trace(inv(A0)*A)         | -KL(A)
%         - trace(WW)                       | p(W)
%
% Let us suppose that WW, C and S are diagonal (which should be the case at
% the optimum). Then the transform we want to optimise is just a diagonal
% scaling matrix, and all individual parameters are independent.
% Let us note for a given component i in [1,K]
%   s         : S(i)
%   c         : C(i)
%   w         : WW(i)
%   a, (N+n0) : A(i)
%   a0, n0    : A0(i)
%   exp(q)    : Q(i) -> Scaling (insure positiveness)
%
% The transform W*inv(Q)*Q*Z implies
%   c <- c * exp(2*q)
%   s <- s * exp(2*q)
%   w <- w * exp(-2*q)
%   a <- update(n0, a0, n0+N, c*exp(2*q))
%
% The update equation for Az writes
%   a = (n0+N) / (n0/a0 + c*exp(2*q))
%
% We thus want to maximise, for each i
%   2*Li = - (n0+N)*c*exp(2*q) / (n0/a0 + c*exp(2*q))
%          - s * exp(2*q)
%          + (n0+1) * (n0+N) / (n0/a0 + c*exp(2*q))
%          - (n0+N) * (n0+N) / (n0 + a0*c*exp(2*q))
%          - w*exp(-2*q)
% 

function solve_symbolic

% compute gradient and hessian
% ----------------------------
% specify positive constraints
% syms sl zl wl n0l Nl a0l q
% n0 = exp(n0l);
% N  = exp(Nl);
% s  = exp(sl);
% z  = exp(zl);
% w  = exp(wl);
% a0 = exp(a0l);
syms s z w n0 N a0 q

c = s + z;
% L = - (n0+N)*c*exp(2*q) / (n0/a0 + c*exp(2*q)) ...
%     - s*exp(2*q) ...
%     + (n0+1) * (n0+N) / (n0/a0 + c*exp(2*q)) ...
%     - (n0+N) * (n0+N) / (n0 + a0*c*exp(2*q)) ...
%     - w*exp(-2*q);
L = - (n0+N)*c*exp(2*q) / (n0/a0 + c*exp(2*q)) ...
    - s*exp(2*q) ...
    + (n0+N) / (n0/a0 + c*exp(2*q)) ...
    - w*exp(-2*q);
L = simplify(L);
G = diff(L, q);
H = diff(G, q);
dq = simplify(-H\G);

% substitute matrix for array operations
% --------------------------------------
dir = char(dq);
dir = strrep(dir, '*', '.*');
dir = strrep(dir, '/', './');
dir = strrep(dir, '^', '.^');
dir

lb = char(L);
lb = strrep(lb, '*', '.*');
lb = strrep(lb, '/', './');
lb = strrep(lb, '^', '.^');
lb

g = char(G);
g = strrep(g, '*', '.*');
g = strrep(g, '/', './');
g = strrep(g, '^', '.^');
g

h = char(H);
h = strrep(h, '*', '.*');
h = strrep(h, '/', './');
h = strrep(h, '^', '.^');
h

end

%%
% Numeric test

function numeric_test

w  = 1;
s  = 1.9;
% z  = 226.6;
z  = 10000;
a0 = 1;
n0 = 32;
N  = 20;
q  = -15:0.1:15;
%lb = ((N + n0).*(n0 + 1))./(n0./a0 + exp(2.*q).*(s + z)) - w.*exp(-2.*q) - (N + n0).^2./(n0 + a0.*exp(2.*q).*(s + z)) - s.*exp(2.*q) - (exp(2.*q).*(N + n0).*(s + z))./(n0./a0 + exp(2.*q).*(s + z));
lb = (N + n0)./(n0./a0 + exp(2.*q).*(s + z)) - w.*exp(-2.*q) - s.*exp(2.*q) - (exp(2.*q).*(N + n0).*(s + z))./(n0./a0 + exp(2.*q).*(s + z));
g = 2.*w.*exp(-2.*q) - 2.*s.*exp(2.*q) - (2.*exp(2.*q).*(N + n0).*(s + z))./(n0./a0 + exp(2.*q).*(s + z)) - (2.*exp(2.*q).*(N + n0).*(s + z))./(n0./a0 + exp(2.*q).*(s + z)).^2 + (2.*exp(4.*q).*(N + n0).*(s + z).^2)./(n0./a0 + exp(2.*q).*(s + z)).^2;
h = (12.*exp(4.*q).*(N + n0).*(s + z).^2)./(n0./a0 + exp(2.*q).*(s + z)).^2 - 4.*w.*exp(-2.*q) - (4.*exp(2.*q).*(N + n0).*(s + z))./(n0./a0 + exp(2.*q).*(s + z)) - (4.*exp(2.*q).*(N + n0).*(s + z))./(n0./a0 + exp(2.*q).*(s + z)).^2 - 4.*s.*exp(2.*q) + (8.*exp(4.*q).*(N + n0).*(s + z).^2)./(n0./a0 + exp(2.*q).*(s + z)).^3 - (8.*exp(6.*q).*(N + n0).*(s + z).^3)./(n0./a0 + exp(2.*q).*(s + z)).^3;
dq = -(2.*s.*exp(2.*q) - 2.*w.*exp(-2.*q) + (2.*exp(2.*q).*(N + n0).*(s + z))./(n0./a0 + exp(2.*q).*(s + z)) + (2.*exp(2.*q).*(N + n0).*(s + z))./(n0./a0 + exp(2.*q).*(s + z)).^2 - (2.*exp(4.*q).*(N + n0).*(s + z).^2)./(n0./a0 + exp(2.*q).*(s + z)).^2)./(4.*s.*exp(2.*q) + 4.*w.*exp(-2.*q) + (4.*exp(2.*q).*(N + n0).*(s + z))./(n0./a0 + exp(2.*q).*(s + z)) + (4.*exp(2.*q).*(N + n0).*(s + z))./(n0./a0 + exp(2.*q).*(s + z)).^2 - (12.*exp(4.*q).*(N + n0).*(s + z).^2)./(n0./a0 + exp(2.*q).*(s + z)).^2 - (8.*exp(4.*q).*(N + n0).*(s + z).^2)./(n0./a0 + exp(2.*q).*(s + z)).^3 + (8.*exp(6.*q).*(N + n0).*(s + z).^3)./(n0./a0 + exp(2.*q).*(s + z)).^3);
[M,I] = max(lb);
qmax = q(I);
fprintf('Max = %6.3g at q = %6.3f\n', M, qmax);
figure
subplot(2,3,1), semilogy(q,lb),        title('objective');
subplot(2,3,2), plot(q,g),         title('gradient');
subplot(2,3,3), semilogy(q,h),         title('hessian');
subplot(2,3,4), plot(q,dq),        title('step');
subplot(2,3,5), plot(q,q+dq),      title('new pos');
subplot(2,3,6), semilogy(q,abs(q+dq-qmax)), title('diff from max');

end