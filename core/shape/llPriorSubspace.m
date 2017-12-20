function ll = llPriorSubspace(W, WW, vs, prm, bnd)
% FORMAT ll = lbSubspace(EW, EWW, (ld / vs, prm, bnd))
% W   - MAP subspace W
% WW  - MAP precision W'*L*W
% vs  - Voxel size
% prm - Parameters of the differential operator
% bnd - Boundary condition (0/1/2/3) [0]
% ld  - LogDet of the operator [recompute]
%
% Log-likelihood of the prior term on the principal subspace
% > ln p(W)

    dim = [size(W) 1 1 1];
    lat = dim(1:3);
    K   = dim(5);
    
    if nargin > 3
        if nargin < 5
            bnd = 0;
        end
        ld = proba('LogDetDiffeo', lat, vs, prm, bnd);
    else
        ld = vs;
    end

    ll = - prod(dim) * log(2*pi) ...
         + K * ld ...
         - trace(numeric(WW));
     ll = ll * 0.5;

end