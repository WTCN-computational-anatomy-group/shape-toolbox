function f = generateObservation(model, opt)

    if isstruct(model.z)
        A   = model.z.A;
        w   = model.pg.w;
        vs  = opt.tpl.vs;
        prm = opt.pg.prm;
        iscat = opt.tpl.cat && checkarray(model.tpl.a);
        if iscat
            a = model.tpl.a;
        else
            a = model.tpl.mu;
        end
    else
        A   = model.Az;
        w   = model.w;
        vs  = opt.vs;
        prm = opt.prm;
        iscat = opt.tpm;
        if iscat && checkarray(model.a)
            a = model.a;
        else
            a = model.mu;
        end
    end

    z = mvnrnd(zeros(1,size(A,1)), inv(A));
    v = reconstructVelocity('latent', z, 'subspace', w);
    clear z
    iphi = exponentiateVelocity(v, 'iphi', 'vs', vs, 'prm', prm);
    clear v
    if iscat
        f = warp(iphi, a);
        f = reconstructProbaTemplate(f);
    else
        f = warp(iphi, a);
    end

end