function lb = lbLatent(dat, model, opt)
% FORMAT ll = lbLatent(dat, model, opt))
%
% Uses model: Az, N, K, ww, zz, Sz
%        dat: z, S z
%        opt: nz0
%
% Part of the lower-bound encompassing latent coordinates terms
% > E[ln p(z)] - E[ln q(z)]

    if opt.nz0 == 0
        lb = - opt.K*opt.N/2*(log(2*pi) - 1);
        lb = lb + model.wpz(1)*opt.N/2*proba('ELogDetWishart', model.Az/(opt.nz0+opt.N), opt.nz0+opt.N);
        lb = lb + model.wpz(2)*opt.N/2*proba('LogDet', model.ww);
    else
        lb = opt.K*opt.N/2;
        lb = lb + opt.N/2*proba('LogDet', ...
            model.wpz(1) * model.Az + model.wpz(2) * model.ww);
    end
    lb = lb - 0.5 * ( model.wpz(1) * trace((model.zz + model.Sz)*model.Az) + ...
                      model.wpz(2) * trace((model.zz + model.Sz)*model.ww) );
    
    tmp = 0;
    for n=1:opt.N
        if any(any(dat(n).Sz ~= 0))
            tmp = tmp + dat(n).z' * (dat(n).Sz \ dat(n).z) ...
                  + proba('LogDet', dat(n).Sz);
        end
    end
    tmp = tmp/2;
    lb = lb+tmp;
end