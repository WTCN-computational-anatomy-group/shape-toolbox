function ll = llPriorSubspace(W, WW, vs, prm, ld)
% FORMAT ll = lbSubspace(EW, EWW, vs, prm, (ld))
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
    
    if nargin < 5
        ld = proba('LogDetDiffeo', lat, vs, prm);
    end

    ll = - prod(dim) * log(2*pi) ...
         + K * ld ...
         - trace(numeric(WW));
     ll = ll * 0.5;

end