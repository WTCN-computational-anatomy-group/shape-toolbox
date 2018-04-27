function model = shape_init_population(model, opt)
% FORMAT dat = shape_init_population(dat, model, opt)
% dat   - Subject-specific data
% model - Model-specific data
% opt   - Options
%
% All inputs are structures that can either be in memory or on disk in the
% form of a mat file. In the latter case, it is read and, if needed,
% written back.
%--------------------------------------------------------------------------
% Initialise population-specific variables.
%
% These parameters can be provided or initialised from scratch.
%  > For file arrays (w, a), this is specified by the [opt.provided]
%    structure.
%  > For precision parameters, they are always initialised from the
%    provided prior value (opt.z.A0, opt.q.A0, opt.v.l0)
%
% They can be considered parameters to optimise or be set fixed
%  > This is specified by the [opt.optimise] structure
%__________________________________________________________________________
% Copyright (C) 2018 Wellcome Centre for Human Neuroimaging

    % =====================================================================
    % Read input from disk (if needed)
    [model, modelpath, opt] = fileToStruct(model, opt);
    
    % =====================================================================
    % Principal subspace
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
        model.pg.ww = computeWLW(model.pg.w, opt.tpl.vs, opt.pg.prm);
    end
    if opt.optimise.pg.w
        model.pg.n = opt.f.N + opt.v.N;
        model      = lbSubspace(model, opt);
    end
    
    % =====================================================================
    % Affine precision
    if opt.f.N
        model.q.A = opt.q.A0;
        if opt.optimise.q.A && opt.q.Mr
            model.q.n = opt.q.n0 + opt.f.N;
            model     = lbAffinePrior(model, opt);
        end
    end
    
    % =====================================================================
    % Latent precision
    model.z.A = opt.z.A0;
    if opt.optimise.z.A
        model.z.n = opt.z.n0 + opt.f.N + opt.v.N;
        model     = lbLatentPrior(model, opt);
    end
   
    % =====================================================================
    % Mixing weight
    model.mixreg.w = opt.mixreg.w0;
    
    % =====================================================================
    % Residual precision
    model.v.l = opt.v.l0;
    if opt.optimise.v.l && opt.v.n0
        model.v.n = model.mixreg.w(1)*(opt.v.N+opt.f.N) + opt.v.n0;
        model     = lbResidualPrior(model, opt);
    end
    
    % =====================================================================
    % Exit
    model = structToFile(model, modelpath); 
end