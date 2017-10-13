function f = generateObservation(model, opt)

    z = mvnrnd(zeros(1,size(model.A,1)), inv(model.A+model.ww));
    v = reconstructVelocity('latent', z, 'subspace', model.w);
    clear z
    iphi = exponentiateVelocity(v, 'iphi', 'vs', opt.vs, 'prm', opt.prm);
    clear v
    f = warp(iphi, model.mu);

end