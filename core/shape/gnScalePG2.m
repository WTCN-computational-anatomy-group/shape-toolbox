function [Q, iQ, q] = gnScalePG2(ww, ezz, n0, N, q0)
% FORMAT [Q, iQ] = gnScalePG2(ww, ezz, n0, N, q0)
% ww  - W'LW
% ezz - Expected value of sum_n z_n * z_n'.
%       > Must have been orthogonalised before.
% n0  - Number of degrees of freedom of the Wishart prior
% N   - Number of observations
% w   - Weight on the WLW part of the z prior [1]
%
% Gauss-Newton optimisation of the scaling factor between LW and Z
    
    K = size(ezz, 1);
    z = diag(ezz);
    w = diag(ww);
    
    if nargin < 5 || isempty(q0)
        q0    = zeros(K,1)-0.5*log(N);
    end
    
    % --- Initial state
    q = min(max(q0,-10),10);
    Q = diag(exp(q));
    iQ = diag(exp(-q));
    E = obj(ww, ezz, Q, iQ, n0, N);
    
    % --- Gauss-Newton optimisation
    for gnit = 1:100
        
        prevE = E;
        
        dq = dir(q, w, z, n0, N);
        
        % --- Line search
        armijo = 1;
        for lsit=1:10
            nq = q + dq/armijo;
            nQ = diag(exp(nq));
            niQ = diag(exp(-nq));
            nE = obj(ww, ezz, nQ, niQ, n0, N);
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

function e = obj(ww, ezz, Q, iQ, n0, N)
% Compute the objective function    

    ezz = Q*ezz*Q;
    ww  = iQ*ww*iQ;
    A   = spm_prob('Wishart', 'up', N, 0, ezz, eye(4), n0);


    e1 = - trace(ww);
    e2 = N * spm_matcomp('LogDet', A);
    e3 = - trace(ezz*A);
    e  = 0.5 * (e1 + e2 + e3);
%     fprintf('%6g %6g %6g %6g %6g\n', e1, e2, e3, e4, e5);

end

function dq = dir(q, w, z, n0, N)
% Compute line search direction
% (Obtained with the symbolic toolbox)

    dq = -(2.*(n0 + N.*z.*exp(2.*q)).*(N.^3.*z.^2.*exp(6.*q) - n0.^2.*w + n0.^2.*z.*exp(4.*q) + N.^2.*n0.*z.*exp(4.*q) - N.^2.*w.*z.^2.*exp(4.*q) + N.*n0.*z.*exp(4.*q) - 2.*N.*n0.*w.*z.*exp(2.*q)))./(4.*n0.^3.*w + 4.*n0.^3.*z.*exp(4.*q) + 4.*N.*n0.^2.*z.*exp(4.*q) + 4.*N.^2.*n0.^2.*z.*exp(4.*q) - 4.*N.*n0.^2.*z.^2.*exp(6.*q) - 4.*N.^2.*n0.*z.^2.*exp(6.*q) + 4.*N.^3.*n0.*z.^2.*exp(6.*q) + 4.*N.^3.*w.*z.^3.*exp(6.*q) + 12.*N.^2.*n0.*w.*z.^2.*exp(4.*q) + 12.*N.*n0.^2.*w.*z.*exp(2.*q));

end

%% ________________________________________________________________________
%
% Symbolic code to obtain the search direction
%  --------------------------------------------
%
% syms q z w n0 N
% 
% e1 = z*exp(2*q)*(n0+N)/(n0+N*z*exp(2*q));
% e2 = w*exp(-2*q);
% e3 = N*log(n0+N*z*exp(2*q));
% e  = simplify(e1 + e2 + e3);
% g  = simplify(diff(e, q));
% h  = simplify(diff(g, q));
% dq = -h\g;
% 
% dqs = char(dq);
% dqs = strrep(dqs, '*', '.*');
% dqs = strrep(dqs, '/', './');
% dqs = strrep(dqs, '\', '.\');
% dqs = strrep(dqs, '^', '.^');