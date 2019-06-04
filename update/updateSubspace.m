function model = updateSubspace(dat, model, opt)
% FORMAT dat = updateSubspace(dat, model, opt)
% model - Model-specific data
% opt   - Options
%
% All inputs are structures that can either be in memory or on disk in the
% form of a mat file. In the latter case, it is read and, if needed,
% written back.
%--------------------------------------------------------------------------
% Update subspace
%__________________________________________________________________________
% Copyright (C) 2018 Wellcome Centre for Human Neuroimaging

    % =====================================================================
    % Read input from disk (if needed)
    [dat, ~, model, modelpath, opt] = fileToStruct(dat, model, opt);
                 
    % =====================================================================
    % Square matrix
    % > M = inv(E[ZZ] + N/(w*lam)*Id)
    M = model.z.S + model.z.zz + model.pg.n * eye(opt.pg.K) / (model.v.l*model.mixreg.w(1));
    M = spm_matcomp('Inv', M);
    
    % =====================================================================
    % Principal subspace
    % > W = V * Z' * M
    %
    % To save disk and memory, I first compute Z' * M (dim = KxN), and 
    % then matrix multiply with V (KxN multiplications).
    %
    % In a "private multicentric" context, we probably would have to first
    % aggregate v * z' in each centre (dim = Ix3xK), and then aggregate
    % those aggregates in the main script.
    P = model.z.Z' * M;
    for k=1:opt.pg.K
        w1 = zeros([opt.tpl.lat 3], 'single');
        for n=1:numel(dat)
            w1 = w1 + numeric(dat(n).v.v) * P(n,k);
        end
        model.pg.w(:,:,:,:,k) = w1;
    end
    clear M P
    model.pg.ww = computeWLW(model.pg.w, opt.tpl.vs, opt.pg.prm);
    
    % =====================================================================
    % Exit
    model = structToFile(model, modelpath);
end