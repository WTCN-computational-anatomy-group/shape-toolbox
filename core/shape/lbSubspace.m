function ll = lbSubspace(EW, EWW, vs, prm, H, reg)
% FORMAT ll = lbSubspace(EW, EWW, vs, prm, (H, reg))
% EW  - Expected subspace E[W]
% EWW - Expected value E[W'*L*W] = I + E[W]'*L*E[W]
% vs  - Voxel size
% prm - Parameters of the differential operator
% H   - Component-specific hessian Hk -ln q(Wk)
% reg - factor of the L'L part of the Hessian
%
% Part of the lower-bound encompassing subspace terms
% > E[ln p(W)] - E[ln q(W)]

    if nargin < 5
        H = [];
    end

    dim = [size(EW) 1 1 1];
    lat = dim(1:3);
    K   = dim(5);

    ll = prod(dim)/2 ...
        + K/2 * proba('LogDetDiffeo', lat, vs, prm) ...
        - 1/2 * trace(numeric(EWW));
    if checkarray(H)
        for k=1:K
            w1 = single(numeric(EW(:,:,:,:,k)));
            h1 = single(numeric(H(:,:,:,:,k)));
            Hw = pointwise3(h1, w1);
            Hw = Hw + reg(k,k)*spm_diffeo('vel2mom', w1, double([vs prm]));
            ll = ll + 0.5 * w1(:)'*Hw(:);
            ll = ll - 0.5 * proba('LogDetDiffeo', lat, vs, double(reg(k,k)*prm));
            if size(h1, 3) == 1
                h1(:,:,:,3) = 1;   % 2D case
            end
            dt = pointwise3(h1, 'd');
            msk = dt > 0;
            ll = ll - 0.5 * sumall(log(dt(msk)));
        end
    end


%     ll = K*(prod(dim)-1)/2;

end