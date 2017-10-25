function lb = lbAffine(dat, model, opt)
% FORMAT ll = lbAffine(dat, model, opt))
%
% Uses model: Aq, nq0, N, qq, Sq
%        dat: q, Sq 
%
% Part of the lower-bound encompassing latent coordinates terms
% > E[ln p(z)] - E[ln q(z)]

    M = size(model.Aq, 1);

    lb = - 0.5*M*opt.N ...
         - 0.5*opt.N*proba('ELogDetWishart', model.Aq/(opt.nq0+opt.N), opt.nq0+opt.N) ...
         - 0.5*trace((model.qq+model.Sq)*model.Aq);
    
    tmp = 0;
    for n=1:opt.N
        if any(any(dat(n).Sq ~= 0))
            tmp = tmp + dat(n).q' * (dat(n).Sq \ dat(n).q) ...
                  + proba('LogDet', dat(n).Sq);
        end
    end
    tmp = 0.5*tmp;
    lb = lb+tmp;
end