function dat = updateResidualLB(dat, model, opt)
% FORMAT dat = updateResidualLB(dat, model, opt)
% dat   - Subject-specific data
% model - Model-specific data
% opt   - Options
%
% All inputs are structures that can either be in memory or on disk in the
% form of a mat file. In the latter case, it is read and, if needed,
% written back.
%--------------------------------------------------------------------------
% Update residual part for lower bound computation
%__________________________________________________________________________
% Copyright (C) 2018 Wellcome Centre for Human Neuroimaging

    % =====================================================================
    % Read input from disk (if needed)
    [dat, datpath, model, ~, opt] = fileToStruct(dat, model, opt);
    
    % The posterior is a multivariate Gaussian with
    % covariance: S = w * lam * W'LW + A
    % mean:       z = w * lam * S * W'Lv
    
    % =====================================================================
    % Reconstruct "mean" velocity
    wz = reconstructVelocity(...
        'latent',   dat.z.z, ...             % Latent coordinates
        'subspace', model.pg.w, ...          % Principal subspace
        'par',      opt.par.within_subject); % Parallelise stuff? (usually no)
    
    
    % =====================================================================
    % Compute residual
    r = numeric(dat.v.v) - wz;
    clear wz
    
    % =====================================================================
    % Compute momentum
    m = spm_diffeo('vel2mom', r, double([opt.tpl.vs opt.pg.prm]));
    
    % =====================================================================
    % Compute residual part
    dat.v.lb.reg = r(:)' * m(:);
    clear r m
                             
    % =====================================================================
    % Exit
    dat = structToFile(dat, datpath);
end