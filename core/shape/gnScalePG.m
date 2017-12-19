function [Q, iQ] = gnScalePG(WW, EZZ, n0, N)
% FORMAT [Q, iQ] = gnScalePG(ww, ezz, n0, N)
% ww  - W'LW
%       > Must have been orthogonalised before.
% ezz - Expected value of sum_n z_n * z_n'.
%       > Must have been orthogonalised before.
% n0  - Number of degrees of freedom of the Wishart prior
% N   - Number of observations
%
% Gauss-Newton optimisation of the scaling factor between LW and Z
    
    K = size(EZZ, 1);

    q    = zeros(K,1)-0.5*log(N);
    q    = min(max(q,-10),10);  % Heuristic to avoid bad starting estimate
    Q    = diag(exp(q));
    A    = spm_prob('Wishart', 'up', N, 0, Q*EZZ*Q, eye(K), n0); % suffstat update
    E    = 0.5*(trace(Q*EZZ*Q*A) + trace(WW/(Q*Q)));
    %fprintf('\n%d %g %g %g\n', 0, 0.5*trace(Q*ZZ1*Q*A), 0.5*trace(WW1*inv(Q*Q)), E)

    for iter=1:100
        A   = spm_prob('Wishart', 'up', N, 0, Q*EZZ*Q, eye(K), n0);
        oE0 = E;

        for subit=1:10
            R  = A.*EZZ'+A'.*EZZ;
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
            E  = 0.5*(trace(Q*EZZ*Q*A) + trace(WW/(Q*Q)));
           %fprintf('%d %g %g %g\n', iter, 0.5*trace(Q*ZZ1*Q*A), 0.5*trace(WW1*inv(Q*Q)), E)
            if (oE-E)/E < 1e-8, break; end
        end
        if abs(oE0-E)/E < 1e-7, break; end
    end
    iQ = inv(Q);