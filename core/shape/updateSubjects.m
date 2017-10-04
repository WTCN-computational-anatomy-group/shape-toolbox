function dat = updateSubjects(model, dat, opt, clean)
% FORMAT dat = updateSubjects(model, dat, opt, ('clean'))
%
% model - Model structure
% dat   - N-structure containing individual images and associated data
% opt   - Option structure 
% clean - If 'clean', clear/remove all non necessary arrays [false]
%
% Update data structure based on the current model parameters (w, ww) and 
% latent variables (z).

    % --- Default parameter
    if nargin < 4
        clean  = '';
    end
    clean = strcmpi(clean, 'clean');
 
    % --- Detect loop/parallelisation scheme
    if strcmpi(opt.loop, 'subject')
        opt.loop = '';
        par      = opt.par;
    else
        par = 0;
    end
    
    % --- Update each subject
    parfor (j=1:numel(dat), double(par))
        dat(j) = updateOneSubject(model, dat(j), opt, clean);
    end

end

function dat = updateOneSubject(model, dat, opt, clean)

    % Exponentiate initial velocity
    % -----------------------------
    dat.v = reconstructVelocity('latent', dat.z, 'subspace', model.w, ...
        'debug', opt.debug, 'output', dat.v, ...
        'loop', opt.loop, 'par', opt.par);
    dat.iphi = exponentiateVelocity(dat.v, 'iphi', ...
        'itgr', opt.itgr, 'vs', opt.vs, ...
        'prm', opt.prm, 'debug', opt.debug, 'output', dat.iphi);
    latf = [size(dat.f) 1];
    latf= latf(1:3);
    dat.ipsi = reconstructIPsi(eye(4), dat.iphi, ...
        'lat', latf, 'Mf', dat.Mf, 'Mmu', model.Mmu, ...
        'output', dat.ipsi, 'debug', opt.debug);


    % Push image to template space
    % ----------------------------
    latmu = [size(model.mu) 1];
    latmu= latmu(1:3);
    [dat.pf, dat.c] = pushImage(dat.ipsi, dat.f, latmu, ...
        'loop', opt.loop, 'par', opt.par, ...
        'output', {dat.pf, dat.c}, 'debug', opt.debug);

    % Compute log-likelihood
    % ----------------------
    dat.llm = llMatching(opt.model, model.mu, dat.pf, dat.c, ...
        'loop', opt.loop, 'par', opt.par, ...
        'debug', opt.debug);

    dat.llz = llPriorLatent(dat.z, model.ww, 'debug', opt.debug);
    
    
    % Clean
    % -----
    if clean
        toclean = {'iphi', 'ipsi'};
        for i=1:numel(toclean)
            field = toclean{i};
            dat.(field) = rmarray(dat.(field));
        end
    end

end