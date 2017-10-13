function ll = llPriorSubspace(W, WW, vs, prm)
% FORMAT ll = lbSubspace(EW, EWW, vs, prm, (H, reg))
% W   - MAP subspace W
% WW  - MAP precision W'*L*W
% vs  - Voxel size
% prm - Parameters of the differential operator
%
% Log-likelihood of the prior term on the principal subspace
% > ln p(W)

    dim = [size(W) 1 1 1];
    lat = dim(1:3);
    K   = dim(5);

    ll = - prod(dim) * log(2*pi) ...
         + K * proba('LogDetDiffeo', lat, vs, prm) ...
         - trace(numeric(WW));
     ll = ll * 0.5;

end