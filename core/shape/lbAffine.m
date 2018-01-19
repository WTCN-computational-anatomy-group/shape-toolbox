function lb = lbAffine(dat, model, opt)
% FORMAT ll = lbAffine(dat, model, opt))
%
% Uses model: Aq, nq0, N, qq, Sq
%        dat: q, Sq 
%
% Part of the lower-bound encompassing latent coordinates terms
% > E[ln p(z)] - E[ln q(z)]
    
    M = size(model.Aq, 1);
    lb = 0;
    if M == 0
        return
    end

    rind = opt.affine_rind;
    for n=1:opt.N
        if any(any(dat(n).Sq ~= 0))
            lb = lb + spm_matcomp('LogDet', dat(n).Sq(rind,rind));
        end
    end
    
    lb = lb + opt.N * M - trace(model.qq(rind,rind) * model.Aq);
    lb = lb + opt.N * spm_prob('Wishart', 'Elogdet', model.Aq, opt.nq0+opt.N, 'normal');
    lb = lb + opt.N * size(model.Aq,1) - trace(model.Sq(rind,rind) * model.Aq);
    
    lb = 0.5 * lb;
end