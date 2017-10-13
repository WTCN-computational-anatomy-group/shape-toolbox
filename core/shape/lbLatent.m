function lb = lbLatent(dat, model, opt)
% FORMAT ll = lbLatent(dat, model, opt))
%
% Uses model: A, n0, N, K, ww, zz, S
%        dat: z, S 
%
% Part of the lower-bound encompassing latent coordinates terms
% > E[ln p(z)] - E[ln q(z)]


    lb = - opt.K*opt.N/2*(log(2*pi) - 1);
    if model.n0 == 0
        lb = lb + model.wpz(1)*opt.N/2*proba('ELogDetWishart', model.A/(model.n0+opt.N), model.n0+opt.N);
    else
        lb = lb + model.wpz(1)*opt.N/2*proba('LogDet', model.A);
    end
    b = lb + model.wpz(2)*opt.N/2*proba('LogDet', model.ww);
    lb = lb - 0.5 * ( model.wpz(1) * trace((model.zz + model.S)*model.A) + ...
                      model.wpz(2) * trace((model.zz + model.S)*model.ww) );
    
    tmp = 0;
    for n=1:opt.N
        if any(any(dat(n).S ~= 0))
            tmp = tmp + dat(n).z' * (dat(n).S \ dat(n).z) ...
                  + proba('LogDet', dat(n).S);
        end
    end
    tmp = tmp/2;
    lb = lb+tmp;
end