function [Q, iQ, q] = pgra_scale_pg(WW, ZZ, SZ, A0, n0, N, q0)
% FORMAT [Q, iQ] = pgva_scale_pg(ww, zz, Sz, l, A0, n0, N, q0)
% ww  - W'LW
%       > Must have been orthogonalised before.
% zz  - Z*Z'
%       > Must have been orthogonalised before.
% Sz  - cov[Z]
%       > Must have been orthogonalised before.
% n0  - Number of degrees of freedom of the Wishart prior
% N   - Number of observations
%
% Gauss-Newton optimisation of the scaling factor between LW and Z
    
    EZZ = ZZ + SZ;
    K = size(EZZ, 1);

    if nargin < 8 || isempty(q0)
        q0    = zeros(K,1)-0.5*log(N);
    end

    q    = min(max(q0,-10),10);  % Heuristic to avoid bad starting estimate
    Q    = diag(exp(q));
    A    = spm_prob('Wishart', 'up', N, 0, Q*EZZ*Q, A0, n0); % suffstat update
    E    = 0.5*(trace(Q*ZZ*Q*A) + trace(WW/(Q*Q)));
%     fprintf('\n%d %g %g %g %g\n', 0, 0.5*trace(Q*ZZ*Q*A), 0.5*trace(WW*inv(Q*Q)), l*trace(SZ*(Q\(WW/Q))), E)

    for iter=1:100
        A   = spm_prob('Wishart', 'up', N, 0, Q*EZZ*Q, A0, n0);
        oE0 = E;

        for subit=1:10
            R  = A.*ZZ'+A'.*ZZ;
            g1 = Q*R*diag(Q);
            g2 =-2*(Q^2\diag(WW));
            g  = g1+g2;

            H1 = Q*R*Q + diag(g1);
            H2 = 4*(Q^2\WW);
            H  = H1+H2;

            H  = spm_matcomp('LoadDiag', H);
            q  = q - H\g;
            q  = min(max(q,-10),10); % Heuristic to avoid overshoot
            Q  = diag(exp(q));

            oE = E;
            E  = 0.5*(trace(Q*ZZ*Q*A) + trace(WW/(Q*Q)));
%             fprintf('\n%d %g %g %g %g\n', 0, 0.5*trace(Q*ZZ*Q*A), 0.5*trace(WW*inv(Q*Q)), l*trace(SZ*(Q\(WW/Q))), E)
            if (oE-E)/E < 1e-8, break; end
        end
        if abs(oE0-E)/E < 1e-7, break; end
    end
    iQ = inv(Q);
    
% % Code for working out the gradients and Hessians
% q   = sym('q',[3,1],'real');
% Q   = diag(exp(q));
% A   = sym('a',[3,3],'real');
% ZZ1 = sym('x',[3,3],'real');
% y   = sym('y',[3,1],'real');
% WW1 = diag(y);
% S1  = sym('z',[3,3],'real');
% l   = sym('l',[1 1],'real');
% %%
% E   = trace(Q*ZZ1*Q*A) + trace(WW1*inv(Q*Q)) + l*trace(S1*inv(Q)*WW1*inv(Q));
% %%
% pretty(simplify(diff(E,sym('q1')),1000))
% pretty(simplify(diff(diff(E,sym('q1')),sym('q2')),1000))
% pretty(simplify(diff(diff(E,sym('q1')),sym('q1')),1000))
% %%
% g1 =  Q*(A.*ZZ1'+A'.*ZZ1)*diag(Q);
% g2 = -Q^2\diag(WW1)*2;
% g3 = -Q^2\(l*diag(WW1*diag(diag(S1))))*2;
% g  =  g1+g2+g3;
% H1 =  Q*(A'.*ZZ1 + A.*ZZ1')*Q +diag(g1);
% H2 =  4*WW1*Q^(-2);
% H3 =  4*l*(WW1*diag(diag(S1)))*Q^(-2);
% H  =  H1+H2+H3;
% %%
% d1  = simplify(g(1)  -diff(E,sym('q1')),1000)
% d11 = simplify(H(1,1)-diff(diff(E,sym('q1')),sym('q1')),1000)
