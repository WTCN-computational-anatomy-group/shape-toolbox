function [Q, iQ] = gnScalePG(ww, zz, S, N, Q0, iQ0, A0, n0)
% FORMAT [Q, iQ] = gnScalePG(ww, zz, S, N, Q0, iQ0)
% ww  - Prior precision matrix of z (W'*L*W)
% zz  - Second order statistics on E[z] (sum_n E[z_n]E[z_n]')
% S   - Sum of expected covariance of z (sum_n cov[z_n])
% N   - Number of subjects
% Q0  - Initial rotation + scaling matrix obtained by SVD
% iQ0 - Initial inverse rotation + scaling matrix obtained by SVD
%
% Optimise the scaling of the principal subspace, while keeping the
% reconstructed velocities (Wz) untouched. This allows to compensate the
% assumption of conditional independence between the subpsace and latent
% coordinates. Its effect is to "move" some of the magnitude from the
% subspace to the latent coordinates prior.
% This scaling is composed with an initial rotation+scaling obtained by
% SVD. This initial transform is that which maximises ln p(Z | W), i.e. 
% Zhang's prior on Z.
% This function refines the scaling by optimising ln p(W) + ln p(Z | A), 
% i.e., our Wishart prior on Z.

    K   = size(ww,1);
    ww  = numeric(ww);
    zz  = numeric(zz);
    S   = numeric(S);
    Q0  = numeric(Q0);
    iQ0 = numeric(iQ0);

    % Apply initial transform
    zz1  = Q0 * zz * Q0';
    Ezz1 = Q0 * (zz + S) * Q0';
    ww1  = iQ0' * ww * iQ0;
    
    % Initialise scaling
    q    = zeros(K,1) - 0.5 * log(N);
    q    = min(max(q, -10), 10);
    Q    = diag(exp(q));
    A    = ww1;
    E    = 0.5*( trace(Q * zz1 * Q * A) + ...
                 trace(ww1 / (Q*Q) ) );

    for iter=1:100
        A   = precisionZWishart(A0, n0, Ezz1, N);
        oE0 = E;

        for subit=1:10
            R  = A .* zz1' + A' .* zz1;
            g1 = Q * R * diag(Q);
            g2 = -2 * ( Q^2 \ diag(ww1));
            g  = g1 + g2;

            H1 = Q * R * Q + diag(g1);
            H2 = 4 * ( Q^2 \ ww1);
            H  = H1 + H2;

            H = loadDiag(H);
            q  = q - H\g;
            q  = min(max(q, -10), 10);
            Q  = diag(exp(q));

            oE = E;
            E  = 0.5 * ( trace(Q * zz1 * Q * A) + trace(ww1 / (Q*Q)) );
            if (oE-E)/E < 1e-8, break; end
        end
        if abs(oE0-E)/E < 1e-7, break; end
    end
    Q  = Q*Q0;
    iQ = iQ0/Q;

end