function lb = lbLatent(dat, model, opt)
% FORMAT ll = lbLatent(dat, model, opt))
%
% Uses model: Az, N, K, ww, zz, Sz
%        dat: z, S z
%        opt: nz0
%
% Part of the lower-bound encompassing latent coordinates terms
% > -KL(q||p) = E[ln p(z)] - E[ln q(z)]

    A = model.wpz(1) * model.Az + model.wpz(2) * model.ww;

    lb = 0;
    for n=1:opt.N
        if any(any(dat(n).Sz ~= 0))
            lb = lb + proba('LogDet', dat(n).Sz);
        end
    end
    
    lb = lb - trace(model.zz * A);
    if opt.nz0 ~= 0
        lb = lb + opt.N * model.wpz(1) * proba('ELogDetWishart', model.Az/(opt.nz0+opt.N), opt.nz0+opt.N);
        lb = lb + opt.N * model.wpz(2) * proba('LogDet', model.ww);
    else
        lb = lb + opt.N * proba('LogDet', A);
    end
    lb = lb + opt.N * size(A,1) - trace(model.Sz * A);
    
    lb = 0.5 * lb;
end