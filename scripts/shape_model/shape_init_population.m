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
        if opt.ui.verbose
            fprintf('%-27s | \n', 'Init Subspace');
        end
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
    if opt.f.N && opt.q.Mr
        if opt.ui.verbose
            fprintf('%-27s | \n', 'Init Affine prec');
        end
        model.q.A = opt.q.A0;
        if opt.optimise.q.A
            model.q.n = opt.q.n0 + opt.f.N;
            if opt.q.n0 == 0 || ~isfinite(opt.q.n0)
                model.q.LogDetA = spm_matcomp('LogDet', model.q.A);
            else
                model.q.LogDetA = spm_prob('Wishart', ...
                    'ELogDet', model.q.A, model.q.n, 'normal');
            end
            model = lbAffinePrior(model, opt);
        end
    end
    
    % =====================================================================
    % Latent precision
    model.z.A = opt.z.A0;
    if opt.optimise.z.A
        if opt.ui.verbose
            fprintf('%-27s | \n', 'Init Latent prec');
        end
        model.z.n = opt.z.n0 + opt.f.N + opt.v.N;
        if opt.z.n0 == 0 || ~isfinite(opt.z.n0)
            model.z.LogDetA = spm_matcomp('LogDet', model.z.A);
        else
            model.z.LogDetA = spm_prob('Wishart', ...
                'ELogDet', model.z.A, model.z.n, 'normal');
        end
        model = lbLatentPrior(model, opt);
    end
   
    % =====================================================================
    % Mixing weight
        if opt.ui.verbose
            fprintf('%-27s | \n', 'Init mixing weight');
        end
    model.mixreg.w = [opt.mixreg.w0 1-opt.mixreg.w0];
    
    % =====================================================================
    % Residual precision
    model.v.l = opt.v.l0;
    if opt.optimise.v.l && opt.v.n0
        if opt.ui.verbose
            fprintf('%-27s | \n', 'Init Residual prec');
        end
        model.v.n = model.mixreg.w(1)*(opt.v.N+opt.f.N) + opt.v.n0;
        if opt.v.n0 == 0 || ~isfinite(opt.v.n0)
            model.v.LogLambda = log(modelv.l);
        else
            model.v.LogLambda = spm_prob('Gamma', ...
                'ELog', model.v.l, model.v.n, 'normal');
        end
        model = lbResidualPrior(model, opt);
    end
    
    % =====================================================================
    % Template
    if ~opt.tpl.provided && opt.f.N && opt.optimise.tpl.a
        if opt.ui.verbose
            fprintf('%-27s | \n', 'Init Template (0)');
        end
        model.tpl.a = prepareOnDisk(model.tpl.a, [opt.tpl.lat opt.model.nc]);
        for z=1:size(model.tpl.a, 3)
            model.tpl.a(:,:,z,:) = 0;
        end
    end
    
    % =====================================================================
    % Exit
    model = structToFile(model, modelpath); 
end