function [dat, model] = pgva_model_init(dat, model, opt)
% Initialise all variables (that need it) 
% + lower bound stuff
    
    
    % ---------------------------------------------------------------------
    %    Model parameters
    % ---------------------------------------------------------------------
    % These parameters can be provided or initialised from scratch.
    %  > For file arrays (w, a), this is specified by the [opt.provided]
    %    structure.
    %  > For precision parameters, they are always initialised from the
    %    provided prior value (opt.z.A0, opt.q.A0, opt.v.l0)
    % They can be considered parameters to optimise or be set fixed
    %  > This is specified by the [opt.optimise] structre
    % ---------------------------------------------------------------------

    % Principal subspace
    % ------------------
    if ~opt.pg.provided
        if opt.optimise.pg.w
            model.pg.w.dim = [opt.tpl.lat 3 opt.pg.K];
            for k=1:opt.pg.K
                model.pg.w(:,:,:,:,k) = 0;
            end
            model.pg.ww = zeros(opt.pg.K);
        else
            model.pg.w.dim = 0;
            model.pg.ww    = 0;
        end
    else
        model.pg.ww = precisionZ(model.pg.w, opt.tpl.vs, opt.pg.prm);
    end
    if opt.optimise.pg.w
        model.pg.n      = opt.f.N + opt.v.N;
        model.lb.w.val  = llPriorSubspace(model.pg.w, model.pg.n * model.pg.ww, opt.pg.ld + prod(opt.tpl.lat)*3*log(model.pg.n));
        model.lb.w.type = 'll';
        model.lb.w.name = 'Subspace prior';
    end
    
    % Affine precision
    % ----------------
    if opt.f.N
        model.q.A = opt.q.A0;
        if opt.optimise.q.A && opt.q.n0 && opt.q.Mr
            model.q.n = opt.q.n0 + opt.f.N;
            model.lb.Aq.val = -spm_prob('Wishart', 'kl', ...
                                        model.q.A,   model.q.n, ...
                                        opt.q.A0,    opt.q.n0, ...
                                        'normal');
            model.lb.Aq.type = 'kl';
            model.lb.Aq.name = '-KL Affine precision';
        end
    end
    
    % Latent precision
    % ----------------
    model.z.A = opt.z.A0;
    if opt.optimise.z.A && opt.z.n0
        model.z.n = opt.z.n0 + opt.f.N + opt.v.N;
        model.lb.Az.val = -spm_prob('Wishart', 'kl', ...
                                    model.z.A,   model.z.n, ...
                                    opt.z.A0,    opt.z.n0, ...
                                    'normal');
        model.lb.Az.type = 'kl';
        model.lb.Az.name = '-KL Latent precision';
    end
    
    % Mixture prior
    % -------------
    model.mixreg.a = opt.mixreg.a0;
    if opt.optimise.mixreg.a
        model.mixreg.n = opt.mixreg.n0 + 1;
        model.lb.ar.val  = 0;
        model.lb.ar.type = 'kl';
        model.lb.ar.name = '-KL Mixture prior';
    end
    
    % Mixture weight
    % --------------
    model.mixreg.w = [model.mixreg.a (1-model.mixreg.a)];
    if opt.optimise.mixreg.w
        model.lb.wr.val  = 0;
        model.lb.wr.type = 'kl';
        model.lb.wr.name = '-KL Mixture weight';
    end
    
    % Residual precision
    % ------------------
    model.v.l = opt.v.l0;
    if opt.optimise.v.l && opt.v.n0
        model.v.n = model.mixreg.w(1)*(opt.v.N+opt.f.N) + opt.v.n0;
        model.lb.l.val  = -spm_prob('Gamma', 'kl', ...
                                    model.v.l, model.v.n, ...
                                    opt.v.l0,  opt.v.n0, ...
                                    prod(opt.tpl.lat)*3, 'normal');
        model.lb.l.type = 'kl';
        model.lb.l.name = '-KL Residual precision';
    end
    
    % ---------------------------------------------------------------------
    %    Individual parameters
    % ---------------------------------------------------------------------
    % These parameters cannot be provided and must be optimised.
    % ---------------------------------------------------------------------
    
    % Momentum
    % --------
    % It is more efficient to compute the momentum of observed velocity
    % fields once and for all.
    dat = pgva_batch('Momentum', dat, opt);
    
    % Affine coordinates
    % ------------------
    if opt.f.N
        [dat, model]    = pgva_batch('InitAffine', 'zero', dat, model, opt);
        model.lb.q.type = 'kl';
        model.lb.q.name = '-KL Affine';
    end
        
    % Latent coordinates
    % ------------------
    [dat, model] = pgva_batch('InitLatent', opt.z.init, dat, model, opt);
    model.lb.z.type = 'kl';
    model.lb.z.name = '-KL Latent';

    % Velocity
    % --------
    [dat, model] = pgva_batch('InitVelocity', 'zero', dat, model, opt);
    % v1 concerns latent velocity fields
    if opt.f.N
        model.lb.v1.type = 'kl';
        model.lb.v1.name = '-KL Velocity';
    end
    % v2 concerns observed velocity fields
    if opt.v.N
        model.lb.v2.type = 'll';
        model.lb.v2.name = 'Velocity likelihood';
    end
    
    % ---------------------------------------------------------------------
    %    Transforms/Template
    % ---------------------------------------------------------------------
    % These parameters depend on the above parameters
    % ---------------------------------------------------------------------
    
    if opt.f.N

        % Initial push/pull
        % -----------------
        % - We need to push observed images to template space to initialise 
        %   the template and to pull (warp) the template to image space to 
        %   initialise the matching term.
        % - It makes sense always keeping the pushed images/pulled template on
        %   disk as it can also be the case in unified segmentation (images 
        %   are responsibilities in this case).
        % - This step also writes ipsi (= iphi(ixi), inverse complete
        %   transform), if they are supposed to be written on disk. It is
        %   usually the case when match = 'pull', since ipsi is needed to
        %   update all warped templates after template update.
        dat = pgva_batch('InitPush', dat, model, opt);

        % Template
        % --------
        if opt.tpl.cat
            if ~opt.tpl.provided
                model.tpl.a = updateMuML(opt.model, dat, ...
                                         'lat',    opt.tpl.lat,   ...
                                         'par',    opt.split.par, ...
                                         'debug',  opt.ui.debug,  ...
                                         'output', model.tpl.a);
            end
            model.tpl.gmu = templateGrad(model.tpl.a,  ...
                                         opt.tpl.itrp, ...
                                         opt.tpl.bnd,  ...
                                         'debug',  opt.ui.debug, ...
                                         'output', model.tpl.gmu);
            model.tpl.mu = reconstructProbaTemplate(model.tpl.a, ...
                                                    'par',    opt.split.par, ...
                                                    'debug',  opt.ui.debug,  ...
                                                    'output', model.tpl.mu);
        else
            if ~opt.tpl.provided
                model.tpl.mu = updateMuML(opt.model, dat, ...
                                          'lat',    opt.tpl.lat,   ...
                                          'par',    opt.split.par, ...
                                          'debug',  opt.ui.debug,  ...
                                          'output', model.tpl.mu);
            end
            model.tpl.gmu = templateGrad(model.tpl.mu, ...
                                         opt.tpl.itrp, ...
                                         opt.tpl.bnd, ...
                                         'debug',  opt.ui.debug, ...
                                         'output', model.tpl.gmu);
        end

        % Matching term
        % -------------
        if strcmpi(opt.match, 'pull')
            [dat, model] = pgva_batch('InitPull', dat, model, opt);
        else
            [dat, model] = pgva_batch('LB', 'Matching', dat, model, opt);
        end
        if opt.f.N
            model.lb.m.type = 'll';
            model.lb.m.name = 'Matching likelihood';
        end
        
    end % < opt.f.N
        
    % Laplace approximation + lower bound stuff
    % -----------------------------------------
    % To compute the KL divergence of distributions estimated by
    % Gauss-Newton, we need images to be pushed or pulled, which can only
    % be done after velocities were initialised
    [dat, model] = pgva_batch('InitLowerBound', dat, model, opt);
    
    % ---------------------------------------------------------------------
    %    Update lower bound
    % ---------------------------------------------------------------------
    model = updateLowerBound(model);    % Accumulate lower bound parts
end
