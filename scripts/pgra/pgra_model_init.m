function [dat, model] = pgra_model_init(dat, model, opt)
% _________________________________________________________________________
%
%    Initialise all variables that need it (if they are not provided).
%
% -------------------------------------------------------------------------
%
% FORMAT [dat, model] = pgra_model_init(dat, model, opt)
%
%
% - Principal geodesic
%   * model.pg.w:  principal subspace [mx my mz 3 k]
%   * model.pg.ww: principal geodesic prior on latent coordinates [k k]
%
% - Latent coordinates
%   * dat.z.z:    parameters [k]
%   * dat.z.zz:   second order (z * z') [k k]
%   * dat.z.S:    posterior covariance matrix [k k]
%   * dat.z.lb:   lots of lower bound / precision stuff to avoid 
%                 recomputing them when updating model parts.
%   * model.z.A:  Regularisation of the affine part (if needed) [k k]
%   * model.z.q:  Sum of all individual parameters (1st order statistic)
%   * model.z.qq: Sum of all individual parameters (2nd order statistic)
%   * model.z.S:  Sum of all posterior covariance matrices
%   * model.z.Z:  Table of all individual z [k N]
%
% - Rigid-body (or affine)
%   * dat.q.q:    parameters (in the Lie algebra) [Q]
%   * dat.q.qq:   second order (q * q') [Q Q]
%   * dat.q.S:    posterior covariance matrix [Q Q]
%   * dat.q.A:    Rigid/affine matrix (exponentiated from q) [4 4]
%   * dat.q.lb:   lots of lower bound / precision stuff to avoid 
%                 recomputing them when updating model parts.
%   * model.q.A:  Regularisation of the affine part (if needed) [Qr Qr]
%   * model.q.q:  Sum of all individual parameters (1st order statistic)
%   * model.q.qq: Sum of all individual parameters (2nd order statistic)
%   * model.q.S:  Sum of all posterior covariance matrices
%
% - Velocity / Residual field
%   * dat.v.r:     residual field [mx my mz 3]
%   * dat.v.v:     initial velocity (= W*z + r) [mx my mz 3]
%   * dat.v.ipsi:  complete inverse transform (= iphi(ixi)) [nx ny nz 3]
%   * dat.v.lb:    lots of lower bound / precision stuff to avoid 
%                  recomputing them when updating model parts.
%   * model.r.l:   Residual precision magnitude (lambda)
%   * model.r.tr:  Sum of all Tr(Sr\L)
%   * model.r.reg: Sum of all r'*L*r
% _________________________________________________________________________
    
    % ---------------------------------------------------------------------
    %    Model parameters
    % ---------------------------------------------------------------------
    % These parameters can be provided or initialised from scratch.
    %  > For file arrays (w, a), this is specified by the [opt.provided]
    %    structure.
    %  > For precision parameters, they are always initialised from the
    %    provided prior value (opt.z.A0, opt.q.A0, opt.r.l0)
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
        model.lb.w.val  = llPriorSubspace(model.pg.w, model.pg.ww, opt.pg.ld);
        model.lb.w.type = 'll';
        model.lb.w.name = 'Subspace prior';
    end
    
    % Affine precision
    % ----------------
    model.q.A = opt.q.A0;
    if opt.optimise.q.A
        model.lb.Aq.val  = 0;
        model.lb.Aq.type = 'kl';
        model.lb.Aq.name = '-KL Affine precision';
    end
    
    % Latent precision
    % ----------------
    model.z.A = opt.z.A0;
    if opt.optimise.z.A
        model.lb.Az.val  = 0;
        model.lb.Az.type = 'kl';
        model.lb.Az.name = '-KL Latent precision';
    end
    
    % Residual precision
    % ------------------
    model.r.l = opt.r.l0;
    if opt.optimise.r.l
        model.lb.l.val  = 0;
        model.lb.l.type = 'kl';
        model.lb.l.name = '-KL Residual precision';
    end
    
    % ---------------------------------------------------------------------
    %    Individual parameters
    % ---------------------------------------------------------------------
    % These parameters cannot be provided and must be optimised.
    % ---------------------------------------------------------------------
    
    % Affine coordinates
    % ------------------
    if opt.optimise.q.q
        [dat, model]    = pgra_batch('InitAffine', 'zero', dat, model, opt);
        model.lb.q.type = 'kl';
        model.lb.q.name = '-KL Affine';
    end
        
    % Latent coordinates
    % ------------------
    if opt.optimise.z.z
        [dat, model] = pgra_batch('InitLatent', opt.z.init, dat, model, opt);
        model.lb.z.type = 'kl';
        model.lb.z.name = '-KL Latent';
    end

    % Velocity
    % --------
    if opt.optimise.r.r
        [dat, model] = pgra_batch('InitResidual', 'zero', dat, model, opt);
        model.lb.r.type = 'kl';
        model.lb.r.name = '-KL Residual';
    end
    
    
    % ---------------------------------------------------------------------
    %    Transforms/Template
    % ---------------------------------------------------------------------
    % These parameters depend on the above parameters
    % ---------------------------------------------------------------------
    
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
    dat = pgra_batch('InitPush', dat, model, opt);
    
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
        [dat, model] = pgra_batch('InitPull', dat, model, opt);
    else
        [dat, model] = pgra_batch('LB', 'Matching', dat, model, opt);
    end
    model.lb.m.type = 'll';
    model.lb.m.name = 'Matching likelihood';
    
    % Laplace approximation
    % ---------------------
    % To compute the KL divergence of distributions estimated by
    % Gauss-Newton, we need images to be pushed or pulled, which can only
    % be done after velocities were initialised
    [dat, model] = pgra_batch('InitLaplace', dat, model, opt);
    
    % ---------------------------------------------------------------------
    %    Update lower bound
    % ---------------------------------------------------------------------
    model = updateLowerBound(model);  % Accumulate lower bound parts
end
